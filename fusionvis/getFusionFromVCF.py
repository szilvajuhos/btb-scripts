import json
import requests, sys
from intervaltree import Interval, IntervalTree
import click
import svgwrite
from svgwrite import cm, mm, rgb
import re

server = "https://rest.ensembl.org"


class SnpEffParser:
    """
        We are processing an snpEff annotation entry (not a VCF line, but the ANN= stuff)
        and returning with a list of dictionaries
    """

    def __init__(self):
        self._ann_list = list()  # list of annotations
        self._ann_keys = list()

    @property
    def ann_keys(self):
        return self._ann_keys

    @property
    def ann_list(self):
        return self._ann_list

    def add_ann_list(self, ann: dict):
        self._ann_list.append(ann)

    def add_ann_keys(self, ann_line: str):
        """
        Parsing the ##INFO=<ID=ANN line
        """
        entry_list = ann_line.split("'")[1]
        entries = entry_list.split("|")
        for item in entries:
            self._ann_keys.append(item.strip())

    def parse_se_ann(self, line, regexp):
        """
        Annotation entries are delimited by commas, process each entry, and add to a list if contains gene_fusion
        """
        ann_dict = {}
        entries = line.split("ANN=")[1]
        entries = entries.split(",")
        for tr_ann in entries:  # get the transcript items
            if regexp.match(tr_ann):
                # chop up values into a list
                ann_list = tr_ann.split("|")
                # add values into a dict
                for i in range(0, len(self.ann_keys)):
                    ann_dict[self.ann_keys[i]] = ann_list[i]
                # add dict to collection (list)
                self._ann_list.append(ann_dict)


class ExonCoords:
    def __init__(self, chromosome, strand, breakpoint, gene_name, exons: IntervalTree):
        self.chromosome = chromosome
        self.strand = strand
        self.breakpoint = breakpoint
        self.gene_name = gene_name
        self.exons = IntervalTree(exons)

    @classmethod
    def fromTuple(cls, a_tuple):
        return cls(a_tuple[0], a_tuple[1], a_tuple[2], a_tuple[3], a_tuple[4])

    @classmethod
    def copy_without_exons(cls, exc):
        return cls(exc.chromosome, exc.strand, exc.breakpoint, exc.gene_name, IntervalTree())

    @classmethod
    def empty(cls):
        return cls("", 0, -1, "", IntervalTree())

    def print_properties(self):
        print("#########################################")
        print("coordinates :", self.chromosome + ":" + str(self.exons.begin()) + "-" + str(self.exons.end()))
        print("gene        :", self.gene_name)
        print("strand      :", self._strand)
        print("breakpoint  :", self._breakpoint)
        print("exons       :", self._exons)
        print("#########################################")

    def print_as_bed(self):
        chromosome = self.chromosome
        for ex in sorted(self.exons):
            print(chromosome + "\t" + str(ex.begin) + "\t" + str(ex.end))

    @property
    def gene_name(self):
        return self._gene_name

    @gene_name.setter
    def gene_name(self, value):
        self._gene_name = value

    @property
    def chromosome(self):
        return self._chromosome

    @chromosome.setter
    def chromosome(self, value):
        self._chromosome = value

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, value):
        self._strand = value

    @property
    def breakpoint(self):  # int
        return self._breakpoint

    @breakpoint.setter
    def breakpoint(self, value):
        self._breakpoint = value

    @property
    def exons(self):  # IntervalTree()
        return self._exons

    @exons.setter
    def exons(self, exons):
        self._exons = exons

    def begin(self):
        return self.exons.begin()


class SV_Maker:
    """
    Joins coordinates of exons for two genes
    """
    # directions from the breakpoint
    DIR_LEFT = True
    DIR_RIGHT = False

    def __init__(self, p5, p3, start, end):
        """
        """
        self.prime5 = ExonCoords(p5.chromosome, p5.strand, 0, p5.gene_name, p5.exons)
        self.prime3 = ExonCoords(p3.chromosome, p3.strand, 0, p3.gene_name, p3.exons)

        # now assign breakpoints to genes:
        # since using the Manta VCF line it is not yet clear
        # which gene contain which endpoint, we have to assign
        # them separately - only if they are meaningful values
        if start is not None and end is not None and start > 0 and end > 0:
            self.assign_breakpoint(start)
            self.assign_breakpoint(end)

    def assign_breakpoint(self, bp):
        for gene in [self.prime5, self.prime3]:
            gene_coords = gene.exons
            gene_ends = IntervalTree()
            gene_ends.add(Interval(gene_coords.begin(), gene_coords.end()))
            if gene_ends.at(bp):
                gene.breakpoint = bp

    def assign_breakpoint_to_genes(self, bp: tuple):
        chromosome = bp[0]
        for gene in [self.prime5, self.prime3]:
            if chromosome == gene.chromosome:
                gene_extremes = IntervalTree()
                gene_extremes.add(Interval(gene.exons.begin(), gene.exons.end()))
                if gene_extremes.at(bp[1]):
                    gene.breakpoint = bp[1]

    def get_left_part(self, gene: ExonCoords):
        # |-->---!->-->-->--|
        # xxxxxxxx
        # if the breakpoint is in an intron, we have to add a 1-base long interval at the breakpoint
        breakpoint_in_exon = False
        ii = iter(gene.exons)
        while ii and not breakpoint_in_exon:
            breakpoint_in_exon = next(ii).contains_point(gene.breakpoint)
        if not breakpoint_in_exon:
            print("**** intron breakpoint at -> ",gene.breakpoint-1,gene.breakpoint)
            gene.exons.add(Interval(gene.breakpoint-1,gene.breakpoint))
        # when the breakpoint is in an exon, we have to shorten that one
        gene.exons.chop(gene.breakpoint, gene.exons.end())

        return gene.exons

    def get_right_part(self, gene: ExonCoords):
        # |--<---!-<--<--<--|
        #        xxxxxxxxxxxx
        tr_exons = gene.exons
        tr_exons = tr_exons.overlap(gene.breakpoint, gene.exons.end())
        tr_exons.add(Interval(gene.breakpoint, gene.breakpoint + 1))
        return tr_exons

    def print_as_bed(self, chromosome, exs):
        for e in sorted(exs):
            print(chromosome + "\t" + str(e.begin) + "\t" + str(e.end))

    def getFusedPart(self, gene: ExonCoords, direction) -> ExonCoords:
        """
        Breaks the gene coordinates at the breakpoint
        and returs with the left or right truncated part
        """
        tr_exons = None
        if direction == 5:
            if gene.strand > 0:
                tr_exons = self.get_left_part(gene)
            else:
                tr_exons = self.get_right_part(gene)
        else:  # assuming prime 3 - the other way around
            if gene.strand < 0:
                tr_exons = self.get_right_part(gene)
            else:
                tr_exons = self.get_left_part(gene)
        return ExonCoords(gene.chromosome, gene.strand, gene.breakpoint, gene.gene_name, tr_exons)

    def fuse_tandem_genes(self):
        # dealing with tandem repeats:
        # first we have to get parts by strand
        # - strand means we want to have the
        #       right part from the breakpoint for 5'
        #       left part from the breakpoint for 3'
        #       join them by starting with the 5' part,
        #       add the 3' part to its left
        # + strand means we want to have the
        #       left part from the breakpoint for 5'
        #       right part from the breakpoint for 3'
        #       join them by starting with the 5' part
        #       add the 3' part to its right
        prime5part = None
        prime3part = None
        if (self.prime5.strand < 0):
            prime5part = ExonCoords(self.prime5.chromosome, self.prime5.strand,
                                    self.prime5.breakpoint, self.prime5.gene_name,
                                    self.get_right_part(self.prime5))
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_left_part(self.prime3))
        else:
            prime5part = ExonCoords(self.prime5.chromosome, self.prime5.strand, self.prime5.breakpoint,
                                    self.prime5.gene_name, self.get_left_part(self.prime5))
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand, self.prime3.breakpoint,
                                    self.prime3.gene_name, self.get_right_part(self.prime3))
        # now move the 3' part to the 5' part
        p5borders = (prime5part.exons.begin(), prime5part.exons.end())
        p3borders = (prime3part.exons.begin(), prime3part.exons.end())
        # |------5------|
        #                   |------3------|
        # |------5------|
        #           |------3------|
        #
        #                   |------5------|
        #           |------3------|
        #                   |------5------|
        # |------3------|
        # shift = (5'start - 3'start)
        # 3'start = 3'start + shift
        shift = 0
        if prime5part.strand > 0:
            shift = prime5part.exons.begin() - prime3part.exons.begin() + 1
        else:
            shift = prime5part.exons.begin() - prime3part.exons.end()
        # we have to shift 3' only
        shifted3p = IntervalTree()
        for iv in prime3part.exons:
            shifted3p.add(Interval(iv.begin + shift, iv.end + shift))
        shifted5p = prime5part.exons
        # and now shift down the stuff to 0 for SVG
        left_shift = (shifted5p | shifted3p).begin()
        # TODO: DRY it out
        based05p = IntervalTree()
        for iv in shifted5p:
            based05p.add(Interval(iv.begin - left_shift, iv.end - left_shift))
        based03p = IntervalTree()
        for iv in shifted3p:
            based03p.add(Interval(iv.begin - left_shift, iv.end - left_shift))
        return (based05p, based03p)

    def fuse_translocations(self, p5dir, p3dir):
        prime5part = None
        prime3part = None
        # TODO: DRY it out
        if p5dir == self.DIR_LEFT:
            prime5part = ExonCoords(self.prime5.chromosome, self.prime5.strand,
                                    self.prime5.breakpoint, self.prime5.gene_name,
                                    self.get_left_part(self.prime5))
        else:
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_right_part(self.prime3))
        if p3dir == self.DIR_LEFT:
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_left_part(self.prime3))
        else:
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_right_part(self.prime3))

        if p5dir == self.DIR_LEFT and p3dir == self.DIR_LEFT:
            # forward antiparallel
            # we have to turn the reverse strand 3' gene backwards
            prime3part = self.turn_backwards(prime3part)
            # and have to stick it to the 5' part
            # we are ignoring chromosomes this time
            # TODO: DRY it out
            based05p = self.shift_left_to(0, prime5part)
            based03p = self.shift_left_to(based05p.end(), prime3part)
            print(based05p)
            print(based03p)
            return (based05p, based03p)

    def shift_left_to(self, new_base: int, exs: ExonCoords):
        rebased = IntervalTree()
        shift = exs.exons.begin() - new_base
        for item in exs.exons:
            rebased.add(Interval(item.begin - shift, item.end - shift))
        return rebased

    def turn_backwards(self, exs: ExonCoords):
        """
        We have coords like:
        |###|----|####|--|#|
        and want to have something like:
        |#|--|####|----|###|
        :param exs:
        ExonCoords that we want to turn backwards
        :return:
        new IntervalTree() with backwards coordinates
        """
        gene_start = exs.exons.begin()
        gene_end = exs.exons.end()
        new_exons = IntervalTree()
        for item in sorted(exs.exons):
            new_end = gene_end - (item.begin - gene_start)
            new_start = gene_end + gene_start - item.end
            # print(exs.chromosome + "\t" + str(new_start) + "\t" + str(new_end) + "\t"+ exs.gene_name)
            new_exons.add(Interval(new_start, new_end))
        return ExonCoords(exs.chromosome, exs.strand, exs.breakpoint, exs.gene_name, new_exons)

    def print_properties(self):
        print("5' gene:")
        print("strand     :", self.prime5.strand)
        print("breakpoint :", self.prime5.chromosome + ":" + str(self.prime5.breakpoint))
        exons = self.prime5.exons
        start = exons.begin()
        end = exons.end()
        print("gene coords:", self.prime5.chromosome + ":" + str(start) + "-" + str(end))
        print("3' gene:")
        print("strand     :", self.prime3.strand)
        print("breakpoint :", self.prime3.chromosome + ":" + str(self.prime3.breakpoint))
        exons = self.prime3.exons
        start = exons.begin()
        end = exons.end()
        print("gene coords:", self.prime3.chromosome + ":" + str(start) + "-" + str(end))
        self.prime5.print_as_bed()
        self.prime3.print_as_bed()


def get_CDS_coords(ENS_ID):
    # look docs at https://rest.ensembl.org/
    print("Looking up " + ENS_ID)
    ext = "/lookup/id/" + ENS_ID + "?expand=1"
    #    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    #    if not r.ok:
    #        r.raise_for_status()
    #        sys.exit()
    with open(ENS_ID + '.json', 'r') as myfile:
        data = myfile.read()
    obj = json.loads(data)
    # now go through the Transcript list
    transcripts = obj['Transcript']
    chromosome = "chr" + str(obj['seq_region_name'])
    exon_intervals = IntervalTree()
    for trs in transcripts:
        # go through each transcript, and store coordinate intervals
        for exon in trs['Exon']:
            start = exon['start']
            end = exon['end']
            exon_intervals.add(Interval(start, end))
    exon_intervals.merge_overlaps()
    print("exons from ENSEMBL JSON:")
    print_exons_as_bed(chromosome, IntervalTree(sorted(exon_intervals.items())), obj['display_name'])
    return (chromosome, obj['strand'], 0, obj['display_name'], IntervalTree(sorted(exon_intervals.items())))

def print_exons_as_bed(chrom,exons,gene_name):
    for item in exons:
        print(chrom + "\t" + str(item.begin) + "\t" + str(item.end) + "\t" + gene_name)


# we are storing pairs of translocations in this dictionary. It is like when we are searching for
# the breakpoint for MantaBND:75600:3:7:0:0:0:1      A       A]CHR6:108561001], it should be
# {"MantaBND:75600:3:7:0:0:0:0":"G]CHR11:65662561]"}
BND_dict = {}


# This is the surrogate for main(): everything happens here

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--vcf', '-v', type=str, help='VCF file to get SVs generated by Manta', required=True)
@click.option('--svg', '-s', type=str, help='The output SVG file name', required=True)
def print_SV(vcf, svg):
    # comment regexp
    comment_re = re.compile("^#.*")
    # further regexps to dig out fusions
    fusion_re = re.compile(".*gene_fusion.*")
    tandem_re = re.compile(".*DUP\:TANDEM.*")
    transloc_re = re.compile(".*MantaBND.*")
    ann_re = re.compile(".*ID=ANN.*")
    # we are dealing with PASS only
    filter_re = re.compile(".*PASS.*")

    # count of pictures (as there can be more)
    pic_count = 0
    # read the VCF file:
    vcf_file = open(vcf, 'r')
    sep = SnpEffParser()
    for line in vcf_file:
        line = line.rstrip()
        # if it is contains the snpEff "ANN" line
        if ann_re.match(line):
            sep.add_ann_keys(line)
            # print("Keys:", sep.ann_keys)
        # it is PASS, has the annotation "gene_fusion" and certainly not a comment
        if filter_re.match(line) and fusion_re.match(line) and not comment_re.match(line):
            sv_call = line.split("\t")
            # process tandem duplications
            if tandem_re.match(line):
                # parse snpEff annotations, and store fusions
                # sep.parse_se_ann(sv_call[7],fusion_re)
                print("######################## processing tandem duplication ####################### ", sv_call[2])
                start = int(sv_call[
                                1])  # column 2 is the SV starting point in the call - just we do not know yet the name of the gene
                snpEff_ann = sv_call[7].split("|")
                # Munching through the "END=140789598;SVTYPE=DUP;SVLEN=1932669;CIPOS=0,1;CIEND=0,1;HOMLEN=1;HOMSEQ=G;SOMATIC;SOMATICSCORE=85;ANN=<DUP:TANDEM>" string to get 140789598
                end = int(snpEff_ann[0].split(";")[0].replace("END=", ""))
                ENS_IDs = snpEff_ann[4].split("&")

                # for forward strand pairs the 5' end is the gene with higher coordinates
                # for reverse strand pairs it is the gene with lower coordinates
                # for tandem duplication fusions the strands should be the same
                genes_to_join = [ExonCoords.fromTuple(get_CDS_coords(ENS_IDs[0])),
                                 ExonCoords.fromTuple(get_CDS_coords(ENS_IDs[1]))]
                # forward strand cases
                if genes_to_join[0].strand > 0 and genes_to_join[1].strand > 0:
                    if genes_to_join[0].begin() < genes_to_join[
                        1].begin():  # we have to swap them as the first is the 3'
                        genes_to_join = [genes_to_join[1], genes_to_join[0]]
                else:  # negative strand
                    if genes_to_join[0].begin() > genes_to_join[1].begin():  # note the relation sign >
                        genes_to_join = [genes_to_join[1], genes_to_join[0]]

                # now we should have the 5' as the first in the list
                prime_5 = genes_to_join[0]
                prime_3 = genes_to_join[1]
                fusion = SV_Maker(prime_5, prime_3, start, end)
                pic_count = makeSVG(fusion.fuse_tandem_genes(), svg, pic_count)
                print("###############################################################################")
            elif transloc_re.match(line):  # translocations
                print("----------------- processing translocation ---------------------- ",
                      sv_call[0], sv_call[1], sv_call[2])
                # get the annotation part
                snpEff_ann = sv_call[7].split("|")
                # look up whether the other end of the translocation is already stored
                mate_ID = snpEff_ann[0].split(";")[1].replace("MATEID=", "")
                print("Searching ", mate_ID)
                if mate_ID in BND_dict.keys():
                    ENS_IDs = snpEff_ann[4].split("&")
                    genes_to_join = [ExonCoords.fromTuple(get_CDS_coords(ENS_IDs[0])),
                                     ExonCoords.fromTuple(get_CDS_coords(ENS_IDs[1]))]
                    # we can have meaningful fusions only for cases like (B is for 'base')
                    # a) genes are parallel (FF or RR), and the chromosome join is B[mate[ - ]mate]B
                    # b) genes are FR, join is B]mate] - B]mate]
                    # c) genes are RF, join is [mate[B - [mate[B
                    # see VCF documentation "5.4 Specifying complex rearrangements with breakends"
                    if genes_to_join[0].strand == genes_to_join[1].strand:
                        print("Parallel strand mates", sv_call[2], sv_call[4],
                              "with mate ID", mate_ID, BND_dict[mate_ID])
                    else:
                        # find out whether it is B]mate]-B]mate] or [mate[B-[mate[B
                        rev_mate_re = re.compile(".*]CHR*.:[0-9].*]")
                        # for B]mate]-B]mate] we want left for both
                        if rev_mate_re.match(sv_call[4]):
                            print("forward antiparallel strand mates", sv_call[2], sv_call[4],
                                  "with mate ID", mate_ID, BND_dict[mate_ID])
                            # the 5' will be the forward gene
                            if genes_to_join[0].strand > 0:
                                (prime_5, prime_3) = (genes_to_join[0], genes_to_join[1])
                            else:
                                (prime_5, prime_3) = (genes_to_join[1], genes_to_join[0])
                            fusion = SV_Maker(prime_5, prime_3, None, None)
                            # have to find out how the breakpoints are assigned
                            # this is for the one in the VCF line
                            fusion.assign_breakpoint_to_genes((sv_call[0], int(sv_call[1])))
                            # this is for the mate
                            breakpoint = extract_breakpoint(sv_call[4])
                            fusion.assign_breakpoint_to_genes(breakpoint)
                            svg_coords = fusion.fuse_translocations(fusion.DIR_LEFT, fusion.DIR_LEFT)
                            pic_count = makeSVG(svg_coords, svg, pic_count)
                            fusion.print_properties()
                        else:
                            # for [mate[B-[mate[B we want right for both
                            print("reverse antiparallel strand mates", sv_call[2], sv_call[4],
                                  "with mate ID", mate_ID, BND_dict[mate_ID])
                else:
                    BND_dict[sv_call[2]] = sv_call[4]


def makeSVG(fex, svg, pic_count):
    w, h = '100%', '100%'
    outfile = str(pic_count) + "_" + svg
    dwg = svgwrite.Drawing(filename=outfile, size=(w, h), debug=True)
    dwg.add(dwg.rect(insert=(0, 0), size=(w, h), fill='white', stroke='black'))
    shapes = dwg.add(dwg.g(id='shapes', fill='red'))
    shapes = shape_intervals(dwg, shapes, fex[0], 'blue')
    shapes = shape_intervals(dwg, shapes, fex[1], 'red')
    shapes.add(dwg.rect(insert=(fex[0].begin() / 100 * mm, 8 * mm),
                        size=(int(fex[0].begin() - fex[0].end()) / 100 * mm, 2 * mm), fill='blue', stroke_width=0))
    shapes.add(dwg.rect(insert=(0, 8 * mm), size=((fex[1].end()) / 100 * mm, 2 * mm), fill='red', stroke='red'))
    dwg.save()
    print("Fusion picture is at", outfile)
    return pic_count + 1

def shape_intervals(dwg, shapes, itv, color):
    height = 2 * cm;
    for iv in itv:
        s = iv.begin / 100
        width = int(abs(iv.end - iv.begin + 100 ) / 100)
        print("SVG:", s*mm, 0, width*mm, height)
        shapes.add(dwg.rect(insert=(s * mm, 0), size=(width * mm, height),
                            fill=color, stroke=color, stroke_width=1))
    return shapes


def extract_breakpoint(a_bp: str):
    # we are getting something like "A]CHR6:108561001]" or "[CHR8:108561001[T" string
    # and we want to return with a ('chr12':123456) tuple
    chrom = ""
    coords = 0
    if a_bp.startswith("["):  # it is like "[CHR8:108561001[T"
        print("bugger [CHR8:108561001[T")
    else:
        mate_list = re.split(':|]', a_bp)  # will be ['A', 'CHR6', '108561001', '']
        chrom = mate_list[1].lower()
        coords = int(mate_list[2])
    return (chrom, coords)


if __name__ == "__main__":
    print_SV()

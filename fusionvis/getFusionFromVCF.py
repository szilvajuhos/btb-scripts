import json
import requests, sys
from intervaltree import Interval, IntervalTree
import click
import svgwrite
from svgwrite import cm, mm, rgb
import re

server = "https://rest.ensembl.org"

class ExonCoords:
    def __init__(self,chromosome, strand, breakpoint, gene_name, exons:IntervalTree):
        self.chromosome = chromosome
        self.strand = strand
        self.breakpoint = breakpoint
        self.gene_name = gene_name
        self.exons = IntervalTree(exons)

    @classmethod 
    def fromTuple(cls, a_tuple):
        return cls( a_tuple[0], a_tuple[1], a_tuple[2], a_tuple[3], a_tuple[4] )

    def print_properties(self):
        print("#########################################")
        print("coordinates :", self.chromosome + ":" + str(self.exons.begin()) +"-"+ str(self.exons.end()))
        print("gene        :", self.gene_name)
        print("strand      :", self._strand)
        print("breakpoint  :", self._breakpoint)
        print("exons       :", self._exons)
        print("#########################################")

    @property
    def gene_name(self):
        return self._gene_name
    @gene_name.setter
    def gene_name(self,value):
        self._gene_name = value

    @property
    def chromosome(self):
        return self._chromosome
    @chromosome.setter
    def chromosome(self,value):
        self._chromosome = value

    @property
    def strand(self):
        return self._strand
    @strand.setter
    def strand(self,value):
        self._strand = value

    @property
    def breakpoint(self):           # int
        return self._breakpoint
    @breakpoint.setter
    def breakpoint(self,value):
        self._breakpoint = value

    @property
    def exons(self):                # IntervalTree()
        return self._exons
    @exons.setter   
    def exons(self,exons):
        self._exons = exons

class SV_Maker:
    """
    Joins coordinates of exons for two genes
    """
    def __init__(self,p5,p3, start, end):
        """
        """
        self.prime5 = ExonCoords(p5.chromosome, p5.strand, 0, p5.gene_name, p5.exons)
        self.prime3 = ExonCoords(p3.chromosome, p3.strand, 0, p3.gene_name, p3.exons)

        # now assign breakpoints to genes:
        # since using the Manta VCF line it is not yet clear 
        # which gene contain which endpoint, we have to assign 
        # them separately
        self.assign_breakpoint(start)
        self.assign_breakpoint(end)

    def assign_breakpoint(self,bp):
        for gene in [self.prime5, self.prime3]:
            gene_coords = gene.exons
            gene_ends = IntervalTree()
            gene_ends.add(Interval(gene_coords.begin(), gene_coords.end()))
            if gene_ends.at(bp):
                gene.breakpoint = bp

    def get_left_part(self,gene:ExonCoords):
            # |-->---!->-->-->--|
            # xxxxxxxx
            tr_exons = gene.exons
            tr_exons = tr_exons.overlap( gene.exons.begin(), gene.breakpoint)
            tr_exons.add(Interval(gene.breakpoint-1,gene.breakpoint))
            return tr_exons

    def get_right_part(self,gene:ExonCoords):
            # |--<---!-<--<--<--|
            #        xxxxxxxxxxxx
            tr_exons = gene.exons
            tr_exons = tr_exons.overlap( gene.breakpoint, gene.exons.end()) 
            tr_exons.add(Interval(gene.breakpoint,gene.breakpoint+1))
            return tr_exons

    def print_as_bed(self,exs):
        chromosome = self.prime5.chromosome
        for e in sorted(exs):
            print(chromosome + "\t" + str(e.begin) + "\t" + str(e.end))
            
    def getFusedPart(self, gene:ExonCoords, prime:int ) -> ExonCoords:
        """
        Breaks the gene coordinates at the breakpoint
        and returs with the left or right truncated part depending on strand
        """
        tr_exons = None
        if prime == 5:
            if gene.strand > 0:
                tr_exons = self.get_left_part(gene)
            else:
                tr_exons = self.get_right_part(gene)
        else:   # assuming prime 3 - the other way around
            if gene.strand < 0:
                tr_exons = self.get_right_part(gene)
            else:
                tr_exons = self.get_left_part(gene)
        return ExonCoords(gene.chromosome, gene.strand, gene.breakpoint, gene.gene_name, tr_exons)

    def fuse_genes(self):
        # first we have to get parts by strand
        # - strand means we want to have the left part from the breakpoint, 
        # + strand means we want to have the right part
        prime5part = self.getFusedPart(self.prime5,5)
        prime3part = self.getFusedPart(self.prime3,3)
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
            shifted3p.add(Interval(iv.begin+shift, iv.end+shift))
        shifted5p = prime5part.exons
        # and now shift down the stuff to 0 for SVG
        left_shift = (shifted5p|shifted3p).begin()
        # TODO: DRY it out
        based05p = IntervalTree()
        for iv in shifted5p:
            based05p.add(Interval(iv.begin-left_shift, iv.end-left_shift))
        based03p = IntervalTree()
        for iv in shifted3p:
            based03p.add(Interval(iv.begin-left_shift, iv.end-left_shift))
        return (based05p,based03p)

    def print_properties(self):
        print("5' gene:")
        print("strand     :", self.prime5.strand)
        print("breakpoint :", self.prime5.chromosome + ":" + str(self.prime5.breakpoint) )
        exons = self.prime5.exons
        start = exons.begin()
        end = exons.end()
        print("gene coords:", self.prime5.chromosome + ":" + str(start) + "-" + str(end) )
        print("3' gene:")
        print("strand     :", self.prime3.strand)
        print("breakpoint :", self.prime3.chromosome + ":" + str(self.prime3.breakpoint) )
        exons = self.prime3.exons
        start = exons.begin()
        end = exons.end()
        print("gene coords:", self.prime3.chromosome + ":" + str(start) + "-" + str(end) )



def get_CDS_coords(ENS_ID):
    # look docs at https://rest.ensembl.org/
    print("Looking up " + ENS_ID)
    ext = "/lookup/id/" + ENS_ID + "?expand=1"
#    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
#    if not r.ok:
#        r.raise_for_status()
#        sys.exit()
    with open(ENS_ID + '.json', 'r') as myfile:
        data=myfile.read()
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
            exon_intervals.add(Interval(start,end))
    exon_intervals.merge_overlaps()
    return (chromosome, obj['strand'], 0, obj['display_name'], IntervalTree(sorted(exon_intervals.items())) )


# This is the surrogate for main(): everything happens here

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcf',      '-v', type=str, help='VCF file to get SVs generated by Manta', required=True)
@click.option('--svg',      '-s', type=str, help='The output SVG file name', required=True)

def print_SV(vcf, svg):
    # comment regexp
    comment_re = re.compile("^#.*")
    # further regexps to dig out fusions 
    fusion_re = re.compile(".*gene_fusion.*")
    tandem_re = re.compile(".*DUP\:TANDEM.*")
    transloc_re = re.compile(".*MantaBND.*")
    # we are dealing with PASS only
    filter_re = re.compile(".*PASS.*")

    # count of pictures (as there can be more)
    pic_count = 0
    # read the VCF file:
    vcf_file = open(vcf,'r')
    for line in vcf_file:
        # it is PASS, has the annotation "gene_fusion" and certainly not a comment
        if filter_re.match(line) and fusion_re.match(line) and not comment_re.match(line):
            line = line.rstrip()
            sv_call = line.split("\t")
            # process tandem duplications
            if tandem_re.match(line):
                print("processing tandem duplication",sv_call[2])
                start = int(sv_call[1])  # column 2 is the SV starting point in the call - just we do not know yet the name of the gene
                snpEff_ann = sv_call[7].split("|")
                # Munching through the "END=140789598;SVTYPE=DUP;SVLEN=1932669;CIPOS=0,1;CIEND=0,1;HOMLEN=1;HOMSEQ=G;SOMATIC;SOMATICSCORE=85;ANN=<DUP:TANDEM>" string to get 140789598
                end = int(snpEff_ann[0].split(";")[0].replace("END=",""))
                ENS_IDs = snpEff_ann[4].split("&")
                prime_5 = ExonCoords.fromTuple(get_CDS_coords(ENS_IDs[0]))
                prime_3 = ExonCoords.fromTuple(get_CDS_coords(ENS_IDs[1]))
                fusion = SV_Maker(prime_5,prime_3, start, end)
                pic_count = makeSVG(fusion.fuse_genes(), svg, pic_count)
                #makeSVG((IntervalTree([Interval(2900,3000),Interval(2700,2800),Interval(1500,1600)]),IntervalTree([Interval(0,100),Interval(200,300),Interval(1400,1500)])), svg)
                #makeSVG((IntervalTree([Interval(0,100),Interval(200,300),Interval(1400,1500)]),IntervalTree([Interval(1500,1600),Interval(2700,2800),Interval(2900,3000)])), svg)
            elif transloc_re.match(line): # usual translocations
                print("processing translocation", sv_call[2])

def makeSVG(fex, svg, pic_count):
    w, h = '100%', '100%'
    outfile = str(pic_count)+"_"+svg
    dwg = svgwrite.Drawing(filename=outfile, size=(w, h), debug=True)
    dwg.add(dwg.rect(insert=(0,0), size=(w, h), fill='white', stroke='black'))
    shapes = dwg.add(dwg.g(id='shapes', fill='red'))
    shapes = shape_intervals(dwg, shapes, fex[0],'blue')
    shapes = shape_intervals(dwg, shapes, fex[1],'red')
    shapes.add(dwg.rect(insert=(fex[0].begin()/100*mm,8*mm), size=(int(fex[0].begin()-fex[0].end())/100*mm,2*mm), fill='blue',stroke_width=0) )
    shapes.add(dwg.rect(insert=(0,8*mm), size=((fex[1].end())/100*mm,2*mm), fill='red',stroke='red') )
    dwg.save()
    print("Fusion picture is at",outfile)
    return pic_count + 1

def shape_intervals(dwg, shapes, itv, color):
    height = 2*cm;
    for iv in itv:
        s = iv.begin/100
        width = int((iv.end - iv.begin)/100)
        shapes.add( dwg.rect(insert=(s*mm,0), size=(width*mm, height), fill=color, stroke=color, stroke_width=1) )
    return shapes
    

if __name__ == "__main__":
    print_SV()

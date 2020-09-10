import json
import requests, sys
from intervaltree import Interval, IntervalTree
import click

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

        print("prime 5 breakpoint:", self.prime5.breakpoint)
        print("prime 3 breakpoint:", self.prime3.breakpoint)
    
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
            if gene.strand > 0:
                tr_exons = self.get_right_part(gene)
            else:
                tr_exons = self.get_left_part(gene)
        return ExonCoords(gene.chromosome, gene.strand, gene.breakpoint, gene.gene_name, tr_exons)

    def fuse_genes(self):
        # first we have to get parts by strand
        # - strand means we want to have the left part from the breakpoint, 
        # + strand means we want to have the right part
        print(type(self.prime5))
        prime5part = self.getFusedPart(self.prime5,5)
        prime3part = self.getFusedPart(self.prime3,3)
        print("-------------Truncated parts---------------")
        prime5part.print_properties()
        prime3part.print_properties()
        print("-------------------------------------------")

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

def print_SV(vcf):
    # read the VCF file:
    vcf_file = open(vcf,'r')
    for line in vcf_file:
        line = line.rstrip()
        sv_call = line.split("\t")
        if sv_call[4] == "<DUP:TANDEM>":
            start = int(sv_call[1])  # column 2 is the SV starting point in the call - just we do not know yet the name of the gene
            snpEff_ann = sv_call[7].split("|")
            # Munching through the "END=140789598;SVTYPE=DUP;SVLEN=1932669;CIPOS=0,1;CIEND=0,1;HOMLEN=1;HOMSEQ=G;SOMATIC;SOMATICSCORE=85;ANN=<DUP:TANDEM>" string to get 140789598
            end = int(snpEff_ann[0].split(";")[0].replace("END=",""))
            ENS_IDs = snpEff_ann[4].split("&")
            prime_5 = ExonCoords.fromTuple(get_CDS_coords(ENS_IDs[0]))
            prime_3 = ExonCoords.fromTuple(get_CDS_coords(ENS_IDs[1]))
            fusion = SV_Maker(prime_5,prime_3, start, end)
            fusion.print_properties()
            fusion.fuse_genes()

if __name__ == "__main__":
    print_SV()

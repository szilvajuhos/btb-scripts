import json
import requests, sys
from intervaltree import Interval, IntervalTree

server = "https://rest.ensembl.org"

def get_CDS_coords(ENS_ID):
    # look docs at https://rest.ensembl.org/
    print("Looking up " + ENS_ID)
    ext = "/lookup/id/" + ENS_ID + "?expand=1"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    obj = r.json()
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
    return(chromosome, obj['strand'], sorted(exon_intervals.items()))

# read the VCF file:
vcf = open('manta_dup.vcf','r')
for line in vcf:
    line = line.rstrip()
    sv_call = line.split("\t")
    if sv_call[4] == "<DUP:TANDEM>":
        snpEff_ann = sv_call[7].split("|")
        ENS_genes = snpEff_ann[4].split("&")
        for gene in ENS_genes:
            print( get_CDS_coords(gene))


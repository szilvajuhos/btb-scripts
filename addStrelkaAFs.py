#!/usr/bin/python3.6
"""
To add germline and somatic allele frequencies for Strelka VCFs you have to massage the entries a bit
SNVs and indels are saved to a separate file in Strelka. It is determined automagically whether it is an
SNV or indel file

Inputs: 
  -v  --vcf     the VCF to annotate

Results are written to STDOUT
"""

import click
import re

class addStrelkaAFs:

    def __init__(self):
        self.fileType = "UNDEFINED"
        self.header_line = re.compile("^#")
        self.chrom_line = re.compile("^#CHROM")
        self.content_line = re.compile("^##content")
        self.snv_content = re.compile(".*snv.*")
        self.indel_content = re.compile(".*indel.*")
        
        # constants:
        self.SNVFILETYPE = "SNVS"
        self.INDELFILETYPE = "INDELS"

    def printVCF(self,vcf):
        with open(vcf,'r') as vcffile:
            self.printHeader(vcffile)
            for line in vcffile:
                  if self.fileType == self.SNVFILETYPE:
                      self.addAFToSNVs(line)
                  elif self.fileType == self.INDELFILETYPE:
                      self.addAFToIndels(line)
                  else:
                      raise Exception("Does not look like a Strelka file for me")

    def printHeader(self,vcffile):
        for line in vcffile:
            if self.chrom_line.match(line):
                print("##INFO=<ID=SAF,Number=.,Type=String,Description=\"Somatic allele frequency for the variant.\">")
                print(line,end="")
                return  # we are returning to the main part, leaving the file pointer at the first non-header line
            elif self.header_line.match(line):
                # the content line should let us know whether it is a Strelka SNVs or indels file
                # if we could not find this defined line, will leave the fileType UNDEFINED
                if self.content_line.match(line):       # there is a content line
                    if self.snv_content.match(line):    # claims to be an SNVs file
                        self.fileType =  self.SNVFILETYPE
                    elif self.indel_content.match(line):
                        self.fileType = self.INDELFILETYPE
                # print out all header lines
                print(line,end="")
        
    def getGenotypeArrays(self, line):
        line = line.rstrip()
        cols = line.split("\t")
        (ref, alt, gtformat) = (cols[3]+"U", cols[4]+"U", cols[8])
        # we also have to split the format and the tumor columns
        gtcols = gtformat.split(":")
        tumorcols = cols[10].split(":") # assuming TUMOR is always in column 11
        
        return (cols,ref,alt,gtformat,gtcols,tumorcols)

    def addSAFtoInfo(self,tier1RefCounts, tier1AltCounts,info):
        try:
            somaticAF = tier1AltCounts / (tier1AltCounts + tier1RefCounts)
        except ZeroDivisionError:   # quite unlikely variant when everything is zero, but whatever
            somaticAF = 0.0
        # have to add an item into INFO
        return info + ";SAF=" + str(round(somaticAF,2))

    # Somatic SNVs:
    # refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FORMAT/AU)
    # altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FORMAT/TU)
    # tier1RefCounts = First comma-delimited value from $refCounts
    # tier1AltCounts = First comma-delimited value from $altCounts
    # Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
    def addAFToSNVs(self,line):
        (cols,ref,alt,gtformat,gtcols,tumorcols) = self.getGenotypeArrays(line)
        refCounts = tumorcols[ gtcols.index(ref) ]
        altCounts = tumorcols[ gtcols.index(alt) ]
        tier1RefCounts = int(refCounts.split(",")[0])
        tier1AltCounts = int(altCounts.split(",")[0])
        # cols[7] is the INFO column
        cols[7] = self.addSAFtoInfo(tier1RefCounts,tier1AltCounts,cols[7])
        print("\t".join(cols))

    # Somatic indels:
    # tier1RefCounts = First comma-delimited value from FORMAT/TAR
    # tier1AltCounts = First comma-delimited value from FORMAT/TIR
    # Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
    def addAFToIndels(self,line):
        (cols,ref,alt,gtformat,gtcols,tumorcols) = self.getGenotypeArrays(line)
        tarCounts = tumorcols[ gtcols.index("TAR") ]
        tier1RefCounts = int(tarCounts.split(",")[0])

        tirCounts = tumorcols[ gtcols.index("TIR") ]
        tier1AltCounts = int(tarCounts.split(",")[0])
        
        # cols[7] is the INFO column
        cols[7] = self.addSAFtoInfo(tier1RefCounts,tier1AltCounts,cols[7])
        print("\t".join(cols))

# This is the surrogate for main(): everything happens here

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcf',      '-v', type=str, help='VCF file to annotate', required=True)

def annotateVCF(vcf):
    afAnn = addStrelkaAFs()
    # make a dict with coords and ranks
    afAnn.printVCF(vcf)

if __name__ == "__main__":
  annotateVCF()

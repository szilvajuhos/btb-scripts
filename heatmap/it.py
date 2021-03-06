import csv
import click
from intervaltree import Interval, IntervalTree
import re

class IntervalPrinter:
    def __init__(self,file_type,infile,chrom,faidx,step):
        self.chrom = chrom
        # calculate chromosome length from FASTA index:
        self.chrLength = self._getChrLength(faidx)
        #print("Chromosome "+self.chrom+" length is "+str(self.chrLength))
        self.t = IntervalTree()
        self.step = step

        if(file_type == 'bedgraph'):
            with open(infile, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if row[0] == chrom:
                        start = int(row[1])
                        end = int(row[2])
                        data = float(row[3]) 
                        if start == end:
                            end = start + 1
                        self.t.addi( start, end, data)
        if(file_type == 'cnvs'):
            with open(infile, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                chrom = chrom.replace("chr","")
                for row in reader:
                    if row[0] == chrom:
                        start = int(row[1])
                        end = int(row[2])
                        data = float(row[3]) 
                        if start == end:
                            end = start + 1
                        self.t.addi( start, end, data)
        if(file_type == 'ratio'):
            with open(infile, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                start = 0   # beginning of the chromosome
                chrom = chrom.replace("chr","")
                for row in reader:
                    if row[0] == chrom:
                        end = int(row[1])
                        #data = float(row[2])   # ratio value
                        data = float(row[4])    # copy number
                        if start == end:
                            end = start + 1
                        self.t.addi( start, end, data)
                        # update
                        start = end

    def _getChrLength(self,faidx):
        with open(faidx,'r') as idx:
            reader = csv.reader(idx,delimiter='\t')
            for row in reader:
                if row[0] == str(self.chrom):
                    return int(row[1])

    def printLine(self):
        sex_re = re.compile(".*[XY]")
        line = ""
        value = 2 
        for i in range(0, self.chrLength, self.step):
            # default value is 2 for autosomes and we have to correct for sex chromosomes below
            value = 2
            # get all the values overlapping the current interval
            overlap = self.t.overlap(i, i+self.step)
            if len(overlap) != 0:
                # we can have more than one intervals overlapping the current one
                data = []
                for interval_obj in overlap:
                    data.append(interval_obj.data)
                value = max(data) #* (1 if not sex_re.match(self.chrom) else 2)
            line = line + str(value) + ","
        line = line + str(value) 
        print(line)


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--file-type',    '-t',
              type=click.Choice(['bedgraph', 'cnvs', 'ratio', 'cnp'], case_sensitive=False),
              help='Input file type',
              required=True)
@click.option('--infile',       '-f', type=str, help='Input file', required=True)
@click.option('--chrom',        '-c', type=str, help='The chromosome in the BED file to process', required=True)
@click.option('--faidx',        '-i', type=str, help='FASTA index file to get chromosome lengths', required=True)
@click.option('--step',         '-s', type=int, help='stepsize', required=False, default=1000000)


def makeCSV(file_type,infile,chrom,faidx,step):
    csvMaker = IntervalPrinter(file_type,infile,chrom,faidx,step)
    csvMaker.printLine()

if __name__ == "__main__":
  makeCSV()

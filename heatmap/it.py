import csv
import click
from intervaltree import Interval, IntervalTree


class IntervalPrinter:
    def __init__(self,bed,chrom,faidx,step):
        self.chrom = chrom
        # calculate chromosome length from FASTA index:
        self.chrLength = self._getChrLength(faidx)
        #print("Chromosome "+self.chrom+" length is "+str(self.chrLength))
        self.t = IntervalTree()
        self.step = step
        with open(bed, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if row[0] == chrom:
                    start = int(row[1])
                    end = int(row[2])
                    data = float(row[3]) 
                    if start == end:
                        end = start + 1
                    self.t.addi( start, end, data)

    def _getChrLength(self,faidx):
        with open(faidx,'r') as idx:
            reader = csv.reader(idx,delimiter='\t')
            for row in reader:
                if row[0] == str(self.chrom):
                    return int(row[1])

    def printLine(self):
        line = ""
        for i in range(0, self.chrLength, self.step):
            value = 1       # default value to one
            if len(self.t[i]) != 0:
                data = []
                for interval_obj in self.t[i]:
                    data.append(interval_obj.data)
                value = max(data)
                if value <= 0.0:
                    value = 0.000001
            line = line + str(value) + ","
        line = line + "1"
        print(line)


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--bed',      '-b', type=str, help='BED file to process', required=True)
@click.option('--chrom',    '-c', type=str, help='The chromosome in the BED file to process', required=True)
@click.option('--faidx',    '-i', type=str, help='FASTA index file to get chromosome lengths', required=True)
@click.option('--step',     '-s', type=int, help='stepsize', required=False, default=1000000)


def makeCSV(bed,chrom,faidx,step):
    csvMaker = IntervalPrinter(bed,chrom,faidx,step)
    csvMaker.printLine()

if __name__ == "__main__":
  makeCSV()

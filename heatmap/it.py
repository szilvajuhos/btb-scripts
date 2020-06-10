import csv
import click
from intervaltree import Interval, IntervalTree


class IntervalPrinter:
    def __init__(self,bed,length,step):
        self.t = IntervalTree()
        self.length = length
        self.step = step
        with open(bed, 'r') as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                start = int(row[1])
                end = int(row[2])
                data = float(row[3]) 
                if start == end:
                    end = start + 1
                self.t.addi( start, end, data)

    def printLine(self):
        line = ""
        for i in range(0,self.length,self.step):
            value = 1       # default value to one
            if len(self.t[i]) != 0:
                data = []
                for interval_obj in self.t[i]:
                    data.append(interval_obj.data)
                value = max(data)
            line = line + str(value) + ","
        line = line + "1"
        print(line)


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--bed',      '-b', type=str, help='BED file to process', required=True)
@click.option('--length',   '-l', type=int, help='chromosome length', required=True)
@click.option('--step',     '-s', type=int, help='stepsize', required=False, default=1000000)


def makeCSV(bed,length,step):
    csvMaker = IntervalPrinter(bed,length,step)
    csvMaker.printLine()

if __name__ == "__main__":
  makeCSV()

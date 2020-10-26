import matplotlib
from matplotlib import pyplot as plt
import numpy as np


def get_chromosome_sizes(index_file):
    chrom_sizes = []
    with open(index_file, 'r') as ifh:
        # get only lines that are starting with "chr", and are not longer than 5 chars
        for line in ifh:
            line_list = line.split()
            chromosome = line_list[0]
            if len(chromosome) < 6 and chromosome.startswith('chr') and chromosome != 'chrM':
                chrom_sizes.append(int(int(line_list[1])/1000000))
    return chrom_sizes


chrom_sizes = get_chromosome_sizes("/home/szilva/genome/Homo_sapiens_assembly38.fasta.fai")
print(chrom_sizes)

fig = plt.figure(figsize=(64,36))
spec = fig.add_gridspec(ncols=24, nrows=1, width_ratios=chrom_sizes, height_ratios=[1])

for plot in list(range(0, 22)) + ['X', 'Y']:
    print(" ", plot, end='')
    if type(plot) == int:
        chrPlot = fig.add_subplot(spec[0,plot], title=str(plot+1))
    if plot == 'X':
        chrPlot = fig.add_subplot(spec[0, 22], title='X')
    if plot == 'Y':
        chrPlot = fig.add_subplot(spec[0, 23], title='Y')

    plt.tick_params(axis='x', which='both', bottom=False, top=False,
                    labelbottom=False)  # labels along the bottom edge are off
    plt.tick_params(axis='y', which='both', left=False, right=False,
                    labelleft=False)  # labels along the left edge are off
    # switch off labels
    chrPlot.set_xticklabels([])
    chrPlot.set_yticklabels([])

plt.subplots_adjust(wspace=0.0, hspace=0.0)

plt.show()
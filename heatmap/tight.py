#from scipy.cluster.hierarchy import linkage
#from scipy.cluster.hierarchy import dendrogram
#from scipy.spatial.distance import pdist
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
#from numpy import arange
import click
#import os


class MakeHeat:
    """
    We are expecting a bunch of CSV files names as chr1.csv, chr2.csv, ... chrX.csv, chrY.csv 
    """
    def __init__(self, step_size, index_file):
        self._chrom_sizes = self._get_chromosome_sizes(index_file)
        self.chromList = list(range(1,23))+['X','Y']
        self.chroms = self._readCSVs()
        (self.min,self.max) = self._getExtremes()
        print("Min: "+str(self.min))
        print("Max: "+str(self.max))
        # change font size
        plt.rcParams.update({'font.size': 36})
        fig = plt.figure(figsize=(64, 36))
        spec = fig.add_gridspec(ncols=24, nrows=1,
                                width_ratios=self._chrom_sizes,
                                height_ratios=[1])

        print("creating plot for chromosome ", 1, end='')
        chr1 = fig.add_subplot(spec[0, 0], title='1')
        # set colormap
        cm = matplotlib.cm.bwr

        samples = np.genfromtxt('chr1.csv', delimiter=',', usecols=[0], dtype=str)

        # generate data
        X = self.chroms['chr1']
        plt.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=True) 
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)

        plt.axvline(x=40)
        chr1.set_yticks(np.arange(X.shape[0]))
        chr1.set_yticklabels(samples)
        plt.pcolor(X,cmap='bwr', vmin=self.min, vmax=self.max)

        ########## next CHRs ################

        for plot in list(range(1,22))+['X','Y']:
            if type(plot) == int:
                chrPlot = fig.add_subplot(spec[0,plot], title=str(plot+1))
                print(" ", plot+1, end='')
            if plot == 'X':
                chrPlot = fig.add_subplot(spec[0, 22], title='X')
                print(" ", 'X', end='')
            if plot == 'Y':
                chrPlot = fig.add_subplot(spec[0, 23], title='Y')
                print(" ", 'Y', end='')

            X = self.chroms['chr'+str(plot)]

            #if plot == 8:
            #    print(X)
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
            plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
            plt.pcolor(X,cmap='bwr', vmin=self.min, vmax=self.max)
            # switch off labels
            chrPlot.set_xticklabels([])
            chrPlot.set_yticklabels([])

        plt.colorbar()
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        print(" ready. Saving plot ...")
        plt.savefig("heatmap.png")
        print("Figure printed to heatmap.png")
        #plt.show()

    def _readCSVs(self):
        """
         we are going to have a CSV file for each chromosome with sample names like
        
         Sample1,0,10,2,0,3,0,4,0
         Sample2,0,10,2,0,3,0,4,0
         Sample3,0,10,2,0,3,0,4,0

         We are not going to save the sample names (assuming they are in the same order for each chromosome), and
         saving data in a dictionary with chromosomes as keys
        """
        chroms = {}
        for c in self.chromList:
            cnvData = np.genfromtxt('chr'+str(c)+'.csv',delimiter=',')
            cnvData = np.delete(cnvData,0,1)    # delete first column
            #cnvData = np.add(cnvData,3.00)  # remove 2 to normalize to zero
            #cnvData = np.log(cnvData)
            #cnvData = np.subtract(cnvData,2.00)  # remove 2 to normalize to zero
            chroms['chr'+str(c)] = cnvData
        return chroms
    
    def _getExtremes(self):
        """
        Find out minimal and maximal values for heat map generation
        """
        minR = np.min(self.chroms['chr1'])  # initialize with extremes of chr1
        maxR = np.max(self.chroms['chr1'])
        for key in self.chroms.keys():
            if np.min(self.chroms[key]) < minR:
                minR = np.min(self.chroms[key])
            if np.max(self.chroms[key]) > maxR:
                maxR = np.max(self.chroms[key])
        
        print(minR,maxR)
        if abs(minR) > maxR:
            maxR = abs(minR)
        else:
            minR = -abs(maxR)
        print(minR,maxR)
        return(-2,6)
        #return(minR,maxR)

    def _get_chromosome_sizes(self, index_file):
        chrom_sizes = []
        with open(index_file, 'r') as ifh:
            # get only lines that are starting with "chr", and are not longer than 5 chars
            for line in ifh:
                line_list = line.split()
                chromosome = line_list[0]
                if len(chromosome) < 6 and chromosome.startswith('chr') and chromosome != 'chrM':
                    chrom_sizes.append(int(int(line_list[1])/10000000))
        return chrom_sizes


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--step_size', '-s', type=int, help='Stepsize [250000]', required=False, default=250000)
@click.option('--chr_index', '-i', type=str, help='Chromosomes index file', required=True)
def printCSV(step_size, chr_index):
    heatMap = MakeHeat(step_size, chr_index)

if __name__ == "__main__":
  printCSV()


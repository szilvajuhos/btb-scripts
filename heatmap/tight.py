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
    def __init__(self, step_size, index_file, centromeres):
        self.chromList = list(range(1,23))+['X','Y']
        self._chrom_den = 10000000  # divide chromosome choord with this
        self._step_size = step_size
        self._chrom_sizes = self._get_chromosome_sizes(index_file)
        self._chrom_ratios = self._get_chromosome_ratios()
        self._centromeres = self._get_centromeres(centromeres)
        self.chroms = self._readCSVs()
        (self.min,self.max) = self._getExtremes()
        print("Min: "+str(self.min))
        print("Max: "+str(self.max))
        print(self._centromeres)
        print(self._chrom_sizes)
        # change font size
        plt.rcParams.update({'font.size': 24})
        # create the figure
        fig = plt.figure(figsize=(64, 36))
        # make a grid for chromosomes (each column is a chromosome, its width is a ratio)
        spec = fig.add_gridspec(ncols=24, nrows=1,
                                width_ratios=self._chrom_ratios,
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

        # draw a line at the centromere
        plt.axvline(x=self.centromere_at('chr1'), linestyle=':', linewidth=1)
        chr1.set_yticks(np.arange(X.shape[0]))
        # switch off labels
        chr1.set_xticklabels([])
        chr1.set_yticklabels(samples)
        # color bar
        plt.pcolor(X, cmap='bwr', vmin=self.min, vmax=self.max)

        ########## next CHRs ################

        for plot in list(range(1,22))+['X', 'Y']:
            chrom_id = str(plot+1) if type(plot) == int else plot
            if type(plot) == int:
                chrPlot = fig.add_subplot(spec[0,plot], title=chrom_id)
                print(" ", chrom_id, end='')
            if plot == 'X':
                chrPlot = fig.add_subplot(spec[0, 21], title=chrom_id)
                print(" ", 'X', end='')
            if plot == 'Y':
                chrPlot = fig.add_subplot(spec[0, 22], title=chrom_id)
                print(" ", 'Y', end='')

            X = self.chroms['chr'+chrom_id]
            # draw a line at the centromere
            plt.axvline(x=self.centromere_at('chr'+chrom_id), linestyle=':', linewidth=1)

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
        plt.show()

    def _get_centromeres(self, cf):
        centromere_coords = {}
        with open(cf,'r') as cfh:
            for line in cfh:
                ls = line.split()
                centromere_coords['chr'+ str(ls[0])] = int(int(ls[1])/self._chrom_den)
        return centromere_coords

    def centromere_at(self,chrom):
        """
        Returns with the place of the vertical line for the centromere
        :param chrom: chromosome (1-22, X, Y)
        :return:
        """
        return self.chroms[chrom].shape[1]*self._centromeres[chrom]/self._chrom_sizes[chrom]

    def _get_chromosome_ratios(self):
        ratios = []
        for chrom in self.chromList:
            ratios.append(self._chrom_sizes['chr'+str(chrom)])
        return ratios

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
        chrom_sizes = {}
        with open(index_file, 'r') as ifh:
            # get only lines that are starting with "chr", and are not longer than 5 chars
            for line in ifh:
                line_list = line.split()
                chromosome = line_list[0]
                if len(chromosome) < 6 and chromosome.startswith('chr') and chromosome != 'chrM':
                    chrom_sizes[str(line_list[0])] = int(int(line_list[1])/self._chrom_den)
        return chrom_sizes


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--step_size', '-s', type=int, help='Stepsize [250000]', required=False, default=250000)
@click.option('--chr_index', '-i', type=str, help='Chromosomes index file', required=True)
@click.option('--centromeres', '-c', type=str, help='Middle point of centromeres (chr coord)', required=True)
def printCSV(step_size, chr_index, centromeres):
    heatMap = MakeHeat(step_size, chr_index, centromeres)

if __name__ == "__main__":
  printCSV()


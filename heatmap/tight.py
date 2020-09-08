from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from numpy import arange
import click
import os


class MakeHeat:
    """
    We are expecting a bunch of CSV files names as chr1.csv, chr2.csv, ... chrX.csv, chrY.csv 
    """
    def __init__(self):
        self.chromList = list(range(1,23))+['X','Y']
        self.chroms = self._readCSVs()
        (self.min,self.max) = self._getExtremes()
        print("Min: "+str(self.min))
        print("Max: "+str(self.max))

        fig = plt.figure()
        #fig = plt.figure(figsize=(15,7))

        chr1 = plt.subplot(1, 24, 1,title='chr1')
        # set colormap
        cm = matplotlib.cm.bwr

        samples = np.genfromtxt('chr1.csv',delimiter=',',usecols=(0),dtype=str)

        # generate data
        X = self.chroms['chr1']
        plt.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=True) 
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        chr1.set_yticks(np.arange(X.shape[0]))
        chr1.set_yticklabels(samples)
        plt.pcolor(X,cmap='bwr', vmin=self.min, vmax=self.max)

        ########## next CHRs ################

        for plot in list(range(2,23))+['X','Y']:
            if type(plot) == int:
                chrPlot = plt.subplot(1, 24, plot,title='chr'+str(plot))
            if plot == 'X':
                chrPlot = plt.subplot(1, 24, 23,title='chrX')
            if plot == 'Y':
                chrPlot = plt.subplot(1, 24, 24,title='chrY')
                
            X = self.chroms['chr'+str(plot)]

            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
            plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
            plt.pcolor(X,cmap='bwr', vmin=self.min, vmax=self.max)
            # switch off labels
            chrPlot.set_xticklabels([])
            chrPlot.set_yticklabels([])

        plt.colorbar()
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        plt.savefig("heatmap.png")
        print("Figure printed to heatmap.png")
        plt.show()

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
        return(0,4)
        return(minR,maxR)

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))

def printCSV():
    heatMap = MakeHeat()

if __name__ == "__main__":
  printCSV()


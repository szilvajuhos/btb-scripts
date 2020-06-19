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
    def __init__(self,csv):

        fig = plt.figure()
        #fig = plt.figure(figsize=(15,7))

        chr1 = plt.subplot(1, 24, 1,title='chr1')
        # set colormap
        cm = matplotlib.cm.bwr

        samples = np.genfromtxt('chr1.csv',delimiter=',',usecols=(0),dtype=str)

        # generate data
        X = np.genfromtxt(csv,delimiter=',')
        #
        # we are going to have a CSV file for each chromosome with sample names like
        #
        # Sample1,0,10,2,0,3,0,4,0
        # Sample2,0,10,2,0,3,0,4,0
        # Sample3,0,10,2,0,3,0,4,0
        #
        # So we have to delete the first column like
        X = np.delete(X,0,1)
        X = np.log(X)
        X = X / X.max()

        plt.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=True) 
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        chr1.set_yticks(np.arange(X.shape[0]))
        chr1.set_yticklabels(samples)
        plt.pcolor(X,cmap='bwr')

        ########## next CHR ################

        for plot in range(2,23):
            chrPlot = plt.subplot(1, 24, plot,title='chr'+str(plot))
            X = np.genfromtxt('chr'+str(plot)+'.csv',delimiter=',')
            X = np.delete(X,0,1)
            X = np.log(X)
            X = X / X.max()

            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
            plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
            plt.pcolor(X,cmap='bwr')
            # switch off labels
            chrPlot.set_xticklabels([])
            chrPlot.set_yticklabels([])

        ######### X,Y ##########
        chrX = plt.subplot(1, 24, 23,title='X')
        X = np.genfromtxt('chrX.csv',delimiter=',')
        X = np.delete(X,0,1)
        X = np.log(X)
        X = X / X.max()

        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
        plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
        plt.pcolor(X,cmap='bwr')
        # switch off labels
        chrPlot.set_xticklabels([])
        chrPlot.set_yticklabels([])

        chrY = plt.subplot(1, 24, 24,title='Y')
        X = np.genfromtxt('chrY.csv',delimiter=',')
        X = np.delete(X,0,1)
        X = np.log(X)
        X = X / X.max()

        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
        plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
        plt.pcolor(X,cmap='bwr')
        # switch off labels
        chrPlot.set_xticklabels([])
        chrPlot.set_yticklabels([])



        plt.colorbar()
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        basename = os.path.basename(csv)
        base = basename.split(".")[0]
        plt.savefig(base + ".png")
        print("Figure printed to "+ base + ".png")
        plt.show()

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--csv',      '-c', type=str, help='CSV file to process', required=True)

def printCSV(csv):
    heatMap = MakeHeat(csv)

if __name__ == "__main__":
  printCSV()


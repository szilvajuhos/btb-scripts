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

        fig = plt.figure(figsize=(15,7))

        ax1 = plt.subplot(1, 2, 1)
        # set colormap
        cm = matplotlib.cm.bwr

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

        # get pairwise distance matrix
        pmat = pdist(X, "euclidean")
        # do hiearchical clustering
        linkmat = linkage(pmat)

        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
        plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
        dendrogram(linkmat,orientation='left')
        x0,x1 = ax1.get_xlim()
        y0,y1 = ax1.get_ylim()

        chr1 = plt.subplot(1, 2, 2)
        Y = X

        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
        plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
        plt.pcolor(Y,cmap='bwr')
        # switch off labels
        chr1.set_xticklabels([])
        chr1.set_yticklabels([])
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


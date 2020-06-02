from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from numpy import arange

fig = plt.figure(figsize=(15,7))

ax1 = plt.subplot(1, 3, 1)
cm = matplotlib.cm.Blues
X = np.random.random([6,6])
pmat = pdist(X, "euclidean")
linkmat = linkage(pmat)
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
dendrogram(linkmat,orientation='left')
x0,x1 = ax1.get_xlim()
y0,y1 = ax1.get_ylim()
#ax1.set_aspect((x1-x0)/(y1-y0))

# add a colorbar to the first plot and immediately make it invisible
#cb = plt.colorbar(ax=ax1)
#cb.ax.set_visible(False)

chr1 = plt.subplot(1, 3, 2)
Y = np.random.random([6,6])

plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
#plt.xticks(arange(0.5, 7.5, 1))
plt.pcolor(Y)
# set a colorbar
#cb = plt.colorbar(ax=chr1)
# and disable
#cb.ax.set_visible(False)
# switch off labels
chr1.set_xticklabels([])
chr1.set_yticklabels([])

chr2 = plt.subplot(1, 3, 3)
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # labels along the bottom edge are off
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) # labels along the left edge are off
Y = np.random.random([16,16])
#plt.xticks(arange(0.5, 7.5, 1))
plt.pcolor(Y)
# switch off labels
chr2.set_xticklabels([])
chr2.set_yticklabels([])
# PLOT the colorbar for the last one
plt.colorbar()


plt.subplots_adjust(wspace=0.0, hspace=0.0)

plt.show()



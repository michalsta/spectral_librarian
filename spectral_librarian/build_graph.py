%load_ext autoreload
%autoreload 2
import networkx as nx
from collections import Counter
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
plt.style.use('dark_background')
from heapq import merge

from spectral_librarian.get_data import list_spectra_clusters
from spectral_librarian.non_overlapping_intervals import OpenClosed, OpenOpen
from spectral_librarian.cluster import PeakCluster, peak_distance_clusters
from spectral_librarian.array_ops import arrayfy

dataFolder = Path("/home/matteo/Projects/eubic2020/spectra_clustering/data/eubic-project")
# dataFolder = Path("~/isotopica/data").expanduser()
mzmlFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML"
mgfFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mgf"
maraclustersFile = dataFolder/"MaRaCluster.clusters_p30.tsv"

# couldn't this be done with multiple threads?
spectra_clusters = list_spectra_clusters(mgfFile, maraclustersFile)


# distribution of precursor charges in collections
Counter(c.precursor_charge() for c in spectra_clusters)
spec_cluster = spectra_clusters[0]
# readjust the distance!
initial_distance = 0.1 # start with this
quorum = .9 # choose only peak clusters with more than 90% of peaks
distance = spec_cluster.refine_intracluster_distance(initial_distance, 0.95)
PC = spec_cluster.peak_clusters(quorum, distance)

PC = np.array(list(PC))
max_pc = max(PC, key=lambda pc: len(pc))
max_pc.plot()

L = np.array([min(pc)[0]-distance for pc in PC])
R = np.array([max(pc)[0]+distance for pc in PC])

repr_mz = np.array([pc[len(pc)//2][0] for pc in PC])
open_intervals = OpenOpen(L, R, sorted=True)

spec_cluster.plot(show=False)
for l,r in zip(L,R):
    plt.axvspan(l,r, alpha=.4, color='green')
plt.show()

x, y = list(zip(*Counter(len(p) for p in PC).items()))
plt.scatter(x,y)
plt.show()

# now, start quarying the intervals.
open_intervals[repr_mz/2]

# then, finally, use lightweigh spectrum





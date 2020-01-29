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
from multiprocessing import Pool

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

max_cluster = max(spectra_clusters, key=len)
max_cluster.plot()
distance = max_cluster.refine_intracluster_distance(.1, 0.95)
PC = list(max_cluster.peak_clusters(.9, distance))



L = np.array([min(pc)[0] - distance for pc in PC])
R = np.array([max(pc)[0] + distance for pc in PC])

for l,r in zip(L,R):
    plt.axvspan(l,r, alpha=.4, color='green')
max_cluster.plot(show=False)
plt.show()



x, y = list(zip(*Counter(len(p) for p in PC).items()))
plt.scatter(x,y)
plt.show()


def get_peak_clusters(spec_cluster, 
                      initial_distance=0.1,
                      quorum=.9):
    # choose only peak clusters with more than 95% of peaks
    distance = spec_cluster.refine_intracluster_distance(initial_distance, 0.95)
    # PC = list(spec_cluster.peak_clusters(quorum, distance))
    return len(spec_cluster), distance
    # L = np.array([min(pc)[0] - distance for pc in PC])
    # R = np.array([max(pc)[0] + distance for pc in PC])
    # repr_mz = np.array([pc[len(pc)//2][0] for pc in PC])
    # interval_spectrum = OpenOpen(L, R, sorted=True)
    # return np.array(PC)
    # , repr_mz, interval_spectrum

Counter(len(sc) for sc in spectra_clusters)
p = Pool(7)
distances = p.map(get_peak_clusters, (sc for sc in spectra_clusters if len(sc)>10))
X = pd.DataFrame(distances, columns=('spec_no','distance'))
plt.scatter(X.spec_no, X.distance)
plt.show()

# X.distance.hist(by=X.spec_no, bins=20)
# plt.show()
# it is stupid to esitmate these things like that..
# could be done only on a selection of spectra.
# For them, simply get the PC and get their left-right distance.
# and use that everywhere.

max_pc = max(PCS[100], key=lambda pc: len(pc))
max_pc.plot()




# now, start quarying the intervals.
open_intervals[repr_mz/2]

# then, finally, use lightweigh spectrum





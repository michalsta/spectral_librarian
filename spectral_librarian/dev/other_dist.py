%load_ext autoreload
%autoreload 2
import networkx as nx
from collections import Counter
import numpy as np
from numpy import array as A
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
plt.style.use('dark_background')
from heapq import merge
from multiprocessing import Pool
import heapq

from spectral_librarian.get_data import list_spectra_clusters
from spectral_librarian.non_overlapping_intervals import OpenClosed, OpenOpen
from spectral_librarian.cluster import PeakCluster, peak_distance_clusters
from spectral_librarian.array_ops import arrayfy
from spectral_librarian.get_our_data import spectra_clusters

max_cluster = max(spectra_clusters, key=len)

MZ = A([p.mz for p in max_cluster.peaks()])
dMZ = np.diff(MZ)


for pc in peak_distance_clusters(self.peaks(), initial_distance):
	if len(pc) > len(self)*quorum and pc.peaks_from_different_spectra():

distance = max_cluster.refine_intracluster_distance(.1, 0.95)
PC = list(max_cluster.peak_clusters(.9, distance))
max_pc = max(PC, key=len)



# problem still there
l_pc, r_pc = max_pc.split(233.153)
l_pc.which_spectra()
r_pc.which_spectra()
max_pc.which_spectra()

pure_PC = [pc for pc in PC if pc.peaks_from_different_spectra()]
len(pure_PC)
len(PC)


w = [pc[-1].mz - pc[0].mz for pc in pure_PC]
plt.scatter([pc[0].mz for pc in pure_PC],
			w)
plt.show()

plt.scatter([pc[0].mz for pc in pure_PC],
			[(pc[-1].mz - pc[0].mz)/len(pc) for pc in pure_PC])
plt.show()
w = np.array([(pc[-1].mz - pc[0].mz) for pc in pure_PC])
plt.hist(w[w>0], bins=100)
plt.show()



x,v = arrayfy(Counter([len(pc) for pc in pure_PC if len(pc) > 10]).items())
plt.scatter(x,v)
plt.show()

big_clusters = heapq.nlargest(20, spectra_clusters, key=len)
sc = next(iter(big_clusters))

MZ = np.array([p.mz for p in sc.peaks()])
plt.scatter(MZ[:-1], np.diff(MZ))
plt.show()
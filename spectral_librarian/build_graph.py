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

from spectral_librarian.get_data import list_collections
from spectral_librarian.non_overlapping_intervals import OpenClosed, OpenOpen

dataFolder = Path("/home/matteo/Projects/eubic2020/spectra_clustering/data/eubic-project")
mzmlFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML"
mgfFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mgf"
maraclustersFile = dataFolder/"MaRaCluster.clusters_p30.tsv"

spectra_collections = list_collections(mgfFile, maraclustersFile)
coll = spectra_collections[0]
coll.precursor_charge()

is_sorted = lambda x: np.all(np.diff(x)>0)


Counter(c.precursor_charge() for c in spectra_collections)

def it(s, spec_no):
    """Generate (mz, intensity, spectrum_number)"""
    for mz, i in zip(s.mz, s.intensity):
        yield mz, i, spec_no

def peak_clusters(spectra, thr=0.1):
    """Generate peak clusters (from lightiest to most massive).

    Arguments:
        spectra to cluster.
    """
    all_confs = merge(*(it(s,i) for i,s in enumerate(spectra)))
    c_prev = next(all_confs)
    peak_cluster = [c_prev]
    for c in all_confs:
        if c[0] - thr >= c_prev[0]:
            yield tuple(peak_cluster)
            peak_cluster = []
        peak_cluster.append(c)
        c_prev = c
    yield tuple(peak_cluster) # last peak cluster

thr = 0.1
PC = np.array(list(peak_clusters(coll, thr)))
L = []
R = []
repr_mz = []
for pclust in PC:
    L.append(min(pclust)[0] - thr)
    R.append(max(pclust)[0] + thr)
    # approximate median mass, but not taking into
    # account the intensities.
    repr_mz.append(pclust[len(pclust)//2][0])

repr_mz = np.array(repr_mz)
open_intervals = OpenOpen(L, R, sorted=True)

# now, start quarying the intervals.
open_intervals[repr_mz/2]




# then, finally, use lightweigh spectrum


OpenClosed




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
from spectral_librarian.cluster import PeakCluster, peak_distance_clusters, ppm_dist_clusters
from spectral_librarian.array_ops import arrayfy
from spectral_librarian.get_our_data import spectra_clusters
from spectral_librarian.models import polyfit

sp = max(spectra_clusters, key=len)
MZ = A([p.mz for p in sp.peaks()])
N = A([p.spec_no for p in sp.peaks()])
dMZ = np.diff(MZ)

plt.scatter(MZ[:-1][dMZ>0], (dMZ/MZ[:-1])[dMZ>0]*1e6, s=1,
			c=N[:-1][dMZ>0])
plt.show()



np.quantile(dMZ, (len(MZ)-1-k0)/(len(MZ)-1))
max([len(s.mz) for s in sp])

plt.plot(MZ[:-1], 1.007/MZ[:-1], c='red', alpha=.5)
plt.plot(MZ[:-1], 1.007/(2*MZ[:-1]), c='red', alpha=.5)
plt.scatter(MZ[:-1], dMZ/MZ[:-1], s=1)
plt.show()
# different approach: there should be a jump in diffs of quantiles
N_Q = 1000
dMZ_Q = np.quantile(dMZ[dMZ>0], np.linspace(0,1,N_Q), interpolation='lower')


plt.plot(np.sort(dMZ))
plt.show()
plt.plot(np.arange(len(dMZ_Q)), dMZ_Q)
plt.show()
plt.plot(dMZ_Q, np.arange(len(dMZ_Q)))
plt.show()
plt.plot(dMZ_Q[dMZ_Q<.1], np.arange(len(dMZ_Q))[dMZ_Q<.1])
plt.show()

plt.plot( np.arange(len(dMZ_Q)-1), np.diff(dMZ_Q) )
plt.show()


def spectral_peak_repeats(peak_clusters):
	return Counter(cnt for pc in peak_clusters
				        for cnt in pc.which_spectra().values())

def left_right_intervals(peak_clusters, d=0):
	L = []; R = []
	for pc in peak_clusters:
		L.append(pc[ 0].mz-d)
		R.append(pc[-1].mz+d)
	return np.array(L), np.array(R)

SPR = spectral_peak_repeats
PDC = peak_distance_clusters

L, R = left_right_intervals(peak_distance_clusters(sp.peaks(), .01), .1)

ppm_clust = [cl for cl in ppm_dist_clusters(sp.peaks(), 100.0) if len(cl) > 10]
SPR(ppm_clust)
L, R = left_right_intervals(ppm_clust)
for l,r in zip(L,R):
    plt.axvspan(l,r, alpha=.4, color='green')
sp.plot()

X = pd.DataFrame({'cnt':len(pc), 'mz_L':pc[0].mz, 'mz_R':pc[-1].mz}
				  for pc in ppm_dist_clusters(sp.peaks(), 100.0))


X['d'] = X.mz_R - X.mz_L
X['mz_ave'] = (X.mz_R + X.mz_L)/2.0

plt.scatter(X.mz_ave[X.d>0], X.d[X.d>0])
plt.show()

m1 = polyfit(x=X.mz_ave[X.d>0], y=(X.d/X.mz_ave)[X.d>0], deg=1)
m1.plot()
m1.error_l1()

m2 = polyfit(x=X.mz_ave[X.d>0], y=X.d[X.d>0], deg=10)
m2.error_l1()
m2.plot()


plt.scatter(X.mz_ave[X.d>0], X.d[X.d>0]/X.mz_ave[X.d>0]*1e6, s=1)
plt.show()


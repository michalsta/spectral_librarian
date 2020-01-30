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

from spectral_librarian.non_overlapping_intervals import OpenClosed, OpenOpen
from spectral_librarian.cluster import PeakCluster, peak_distance_clusters, ppm_dist_clusters
from spectral_librarian.array_ops import arrayfy
from spectral_librarian.get_our_data import spectra_clusters
from spectral_librarian.models import polyfit

sp = max(spectra_clusters, key=len)
MZ = A([p.mz for p in sp.peaks()])
N = A([p.spec_no for p in sp.peaks()])
dMZ = np.diff(MZ)

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
LRI = left_right_intervals

L, R = LRI(peak_distance_clusters(sp.peaks(), .01), .1)

ppm_clust = [cl for cl in ppm_dist_clusters(sp.peaks(), 100.0) if len(cl) > 10]
SPR(ppm_clust)
L, R = left_right_intervals(ppm_clust)
for l,r in zip(L,R):
    plt.axvspan(l,r, alpha=.4, color='green')
sp.plot()

abs_clust = [cl for cl in peak_distance_clusters(sp.peaks(), 0.05) if len(cl) > 10]
SPR(peak_distance_clusters(sp.peaks(), 0.05))
SPR(peak_distance_clusters(sp.peaks(), 0.02))
SPR(peak_distance_clusters(sp.peaks(), 0.01))

X = pd.DataFrame({'cnt':len(pc), 
				  'mz_L':pc[0].mz,
				  'mz_R':pc[-1].mz}
				  for pc in ppm_dist_clusters(sp.peaks(), 100.0))

X['d'] = X.mz_R - X.mz_L
X['mz_ave'] = (X.mz_R + X.mz_L)/2.0

plt.scatter(X.mz_ave[X.d>0], X.d[X.d>0])
plt.show()
m1 = polyfit(x=X.mz_ave[X.d>0], y=(X.d/X.mz_ave)[X.d>0], deg=1)
m1.plot()



# def refine_ppm_dist()
# check, if this makes sense
quorum = .2
MZ = []
peak_no = 0
peak_clusters_no = 0
add_clusters = 0
for pc in ppm_dist_clusters(sp.peaks(), 100.0):
    if len(pc) > len(sp)*quorum and pc.peaks_from_different_spectra():
        for p in pc:
            MZ.append(p.mz)
        peak_no += len(pc)
        peak_clusters_no += 1
        add_clusters += max(pc.which_spectra().values()) - 1
MZ = np.array(MZ)
dMZ = np.diff(MZ)

last_peak_diff = (peak_no-1-peak_clusters_no)/(peak_no-1)
quantile_distance = np.quantile(2*dMZ/(MZ[1:]+MZ[:-1])*1e6,
								last_peak_diff)

P = np.linspace( 0,1,10000)
Q = np.quantile(2*dMZ/(MZ[1:]+MZ[:-1])*1e6, P)

plt.plot(Q, P)
plt.scatter(quantile_distance, last_peak_diff)
plt.show()

# what if we modified the k?
MZ = A([p.mz for p in sp.peaks()])
N = A([p.spec_no for p in sp.peaks()])
dMZ = np.diff(MZ)
peak_no = len(MZ)

last_peak_diff = (peak_no-1-peak_clusters_no)

P = np.arange(peak_no)/(peak_no-1)
ppm_d_consec_peaks = 2*dMZ/(MZ[1:]+MZ[:-1])*1e6
Q = np.quantile(ppm_d_consec_peaks, P, interpolation='lower')

np.quantile(ppm_d_consec_peaks, 1-1/(.8*len(sp)), interpolation='lower')

plt.plot(P,Q)
plt.show()

plt.hist(Q[Q<100], bins=1000)
plt.show()

plt.plot(Q,P)
plt.show()

plt.plot(P[:-2], np.diff(Q, 2))
plt.show()

# the triangle method:
plt.plot(P_100, Q_100)
plt.show()

Q_100 = Q[Q <= 100]
P_100 = P[:len(Q_100)]


v1 = Q_100[-1] - Q_100[0]
v2 = P_100[-1] - P_100[0]
a = v2 / v1
b = Q_100[0] - a*P_100[0]
adjusted_P = P_100 - a*Q_100 - b
i = np.argmax(adjusted_P)
Q_100[i]

# v_norm = np.sqrt(v1**2 + v2**2)
# v1, v2 = v1/v_norm, v2/v_norm
# dP_100 = P_100[1] - P_100[0]  
# dQ_100 = np.diff(Q_100)
# norms = np.sqrt(dQ_100**2 + dP_100**2)
# alphas = (dQ_100*v1 + dP_100*v2)/norms
# i = np.argmax(alphas)
plt.plot(Q_100, P_100)
plt.scatter(Q_100[i], P_100[i])
plt.show()


plt.plot(Q_100, P_100)
plt.show()
f = Q_100[-1]
x = P_100[-1]
X = (f/(x**2 + f**2)) * np.array([[x/f, 1], [-f, x]])
Y = np.dot(X, np.array([P,Q]))
plt.plot(Y[0,:], Y[1,:])
plt.show()

i = np.argmin(Y[1,:])

plt.plot(Q_100, P_100)
plt.scatter(Q_100[i], P_100[i])
plt.show()




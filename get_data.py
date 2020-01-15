from pathlib import Path
import pymzml
import json
import numpy as np

import IsoSpecPy

datafolder = Path("~/Projects/eubic2020/spectra_clustering/data/eubic-project").expanduser()
mzml_file = list((datafolder).glob('*.mzML'))[0]

mzml = pymzml.run.Reader(str(mzml_file))
MZML = list(mzml)
MZML = np.array(MZML)

path = datafolder.glob

def get_cluster_indices(path):


idx = {2225, 2278, 2292, 2299, 2306, 2334, 
	2285, 
	2320, 
	2313, 
	2370, 
	2327, 
	2348, 
	2363, 
	2229, 
	2254, 
	2236, 
	2250, 
	2271, 
	2261, 
	2243, 
	2341, 
	2355, 
	2359, 
	2390}

cluster0 = MZML[np.array(list(idx))-1]




spec = next(iter(mzml))
spec.scan_time_in_minutes()
spec.ms_level
spec.selected_precursors

import matplotlib.pyplot as plt

spec = cluster0[1]
plt.stem(spec.mz, spec.i/spec.i.sum())

spec = cluster0[10]
plt.stem(spec.mz, -spec.i/spec.i.sum())

plt.show()

s = IsoSpecPy.IsoDistribution(masses=spec.mz, probs=spec.i)
s.plot(show=True)




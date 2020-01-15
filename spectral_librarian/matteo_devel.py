%load_ext autoreload
%autoreload 2

import IsoSpecPy
from pathlib import Path
import json
import numpy as np

from spectral_librarian.get_data import parse_maracluster_output, get_spectra_array
from spectral_librarian.plot import plot_spectrum
from spectral_librarian.model import SpectralModel

datafolder = Path("~/Projects/eubic2020/spectra_clustering/data/eubic-project").expanduser()

mzml_file = list((datafolder).glob('*.mzML'))[0]
MZML = get_spectra_array(mzml_file)

maracluter_outcome_path = list(datafolder.glob('*.tsv'))[0]
clusters_idx = parse_maracluster_output(maracluter_outcome_path)

clusters = [list(MZML[i-1]) for i in clusters_idx]
spec = clusters[0][1]
# plot_spectrum(spec[0], spec[1])

res = SpectralModel.fromMZPList(clusters[0])

import json

with open('clust0.json', 'w') as f:
	json.dump([tuple(s) for s in clusters[0]], f)

for s in clusters[0]

clusters[0][0]
len(res.peak_clusters)

# spec = next(iter(mzml))
# spec.scan_time_in_minutes()
# spec.ms_level
# spec.selected_precursors

# import matplotlib.pyplot as plt

# spec = cluster0[1]
# plt.stem(spec.mz, spec.i/spec.i.sum())

# spec = cluster0[10]
# plt.stem(spec.mz, -spec.i/spec.i.sum())

# plt.show()

# s = IsoSpecPy.IsoDistribution(masses=spec.mz, probs=spec.i)
# s.plot(show=True)




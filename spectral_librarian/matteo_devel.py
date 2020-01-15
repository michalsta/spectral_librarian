%load_ext autoreload
%autoreload 2

import IsoSpecPy
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt

from spectral_librarian.get_data import parse_maracluster_output, get_spectra_array
from spectral_librarian.plot import plot_spectrum, side_plot_spectra
from spectral_librarian.model import SpectralModel

datafolder = Path("~/Projects/eubic2020/spectra_clustering/data/eubic-project").expanduser()

mzml_file = list((datafolder).glob('*.mzML'))[0]
MZML = get_spectra_array(mzml_file)
maracluter_outcome_path = list(datafolder.glob('*.tsv'))[0]
clusters_idx = parse_maracluster_output(maracluter_outcome_path)


clusters = [[MZML[i-1] for i in I] for I in clusters_idx]
# spec = clusters[0][1]
# plot_spectrum(spec[0], spec[1])

res = SpectralModel.fromMZPList(clusters[0])
res.peak_clusters

res.plot()

peak_clusters = res.peak_clusters
cl = peak_clusters[49]
it = (cl for cl in peak_clusters if len(cl) > 4)
cl = next(it)

with open('clust0.json', 'w') as f:
    json.dump(clusters[0], f)




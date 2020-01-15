%load_ext autoreload
%autoreload 2

import IsoSpecPy
from pathlib import Path
import json
import numpy as np

from get_data import parse_maracluster_output, get_spectra_array
from plot import plot_spectrum

datafolder = Path("~/Projects/eubic2020/spectra_clustering/data/eubic-project").expanduser()

mzml_file = list((datafolder).glob('*.mzML'))[0]
MZML = get_spectra_array(mzml_file)

maracluter_outcome_path = list(datafolder.glob('*.tsv'))[0]
clusters_idx = parse_maracluster_output(maracluter_outcome_path)

clusters = [MZML[i-1] for i in clusters_idx]
clusters[0][1]
spec

from plot import plot_spectrum




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




from pathlib import Path

from spectral_librarian.get_data import list_spectra_clusters

dataFolder = Path("/home/matteo/Projects/eubic2020/spectra_clustering/data/eubic-project")
# dataFolder = Path("~/isotopica/data").expanduser()
mzmlFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML"
mgfFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mgf"
maraclustersFile = dataFolder/"MaRaCluster.clusters_p30.tsv"
spectra_clusters = list_spectra_clusters(mgfFile, maraclustersFile)

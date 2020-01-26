from collections import defaultdict
import numpy as np

# Wouts (most likely) scripts
from spectral_librarian.ms_io import read_spectra, read_clusters, write_spectra, read_psms
from spectral_librarian.cluster import SpectraCluster


def list_spectra_clusters(mgfFile, maraclustersFile):
    """Get a list of spectra clusters made with MaRaCluster.

    Arguments:
        mgfFile (str,Path): path to mgf file with spectra.
        maraclustersFile (str,Path): path to the MaRaCluster clusters.

    Returns:
        list of spectra clusters.
    """
    spectra = {f'{s.filename}:scan:{s.scan}': s
               for s in read_spectra(str(mgfFile))}
    spec2clust = read_clusters(str(maraclustersFile), 'maracluster')
    spectra_clusters = defaultdict(SpectraCluster)
    for sp, cl in spec2clust.items():
        spectra_clusters[cl].append(spectra[sp])
    return list(spectra_clusters.values())# list of Collections of spectra



if __name__ == '__main__':
    import pandas as pd
    from pathlib import Path
    from collections import Counter
    import matplotlib.pyplot as plt

    # input
    dataFolder = Path("/home/matteo/Projects/eubic2020/spectra_clustering/data/eubic-project")
    mzmlFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML"
    mgfFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mgf"
    maraclustersFile = dataFolder/"MaRaCluster.clusters_p30.tsv"

    print('Run it with something higher than Python3.5')
    spectra_clusters = list_spectra_clusters(mgfFile, maraclustersFile)

    # all collections have the same charge!
    collection_with_different_precursor_charges = [sp for sp in spectra_clusters if len(Counter(sp.precursor_charges())) > 1] 
    print(collection_with_different_precursor_charges)
    
    clusters_stats = pd.DataFrame({'size': len(sp),
                                   'prec_mz_mean': sp.precursor_mzs().mean(),
                                   'prec_mz_std': sp.precursor_mzs().std()
                                   } for sp in spectra_clusters)
    plt.scatter(clusters_stats.prec_mz_mean,
                clusters_stats.prec_mz_std)
    plt.xlabel('Mean of the cluster precursor m/z')
    plt.ylabel('Stdev of cluster precursor m/z')
    plt.show()
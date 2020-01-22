from pathlib import Path
from collections import defaultdict, Counter
import spectrum_utils.plot as sup
import matplotlib.pyplot as plt
from itertools import islice
import numpy as np

# Wouts (most likely) scripts
from ms_io import read_spectra, read_clusters, write_spectra, read_psms
import ms_io


class Collection(list):
    """A simple extension of list to store one cluster of spectra."""
    def precursor_charges(self):
        return np.array([s.precursor_charge for s in self])

    def precursor_mzs(self):
        return np.array([s.precursor_mz for s in self])

def list_collections(mgfFile, maraclustersFile):
    """Get a list of collections of spectra clustered by MaRaCluster.

    Arguments:
        mgfFile (str,Path): path to mgf file with spectra.
        maraclustersFile (str,Path): path to the MaRaCluster clusters.

    Returns:
        list of Collection objects.
    """
    spectra = {f'{s.filename}:scan:{s.scan}': s
               for s in ms_io.read_spectra(str(mgfFile))}
    spec2clust = read_clusters(str(maraclustersFile), 'maracluster')
    spectra_collections = defaultdict(Collection)
    for sp, cl in spec2clust.items():
        spectra_collections[cl].append(spectra[sp])
    return list(spectra_collections.values())# list of Collections of spectra


if __name__ == '__main__':
    import pandas as pd

    # input
    dataFolder = Path("/home/matteo/Projects/eubic2020/spectra_clustering/data/eubic-project")
    mzmlFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML"
    mgfFile = dataFolder/"01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mgf"
    maraclustersFile = dataFolder/"MaRaCluster.clusters_p30.tsv"

    print('Run it with something higher than Python3.5')
    spectra_collections = list_collections(mgfFile, maraclustersFile)

    # all collections have the same charge!
    collection_with_different_precursor_charges = [sp for sp in spectra_collections if len(Counter(sp.precursor_charges())) > 1] 
    print(collection_with_different_precursor_charges)
    
    clusters_stats = pd.DataFrame({'size': len(sp),
                                   'prec_mz_mean': sp.precursor_mzs().mean(),
                                   'prec_mz_std': sp.precursor_mzs().std()
                                   } for sp in spectra_collections)
    plt.scatter(clusters_stats.prec_mz_mean,
                clusters_stats.prec_mz_std)
    plt.xlabel('Mean of the cluster precursor m/z')
    plt.ylabel('Stdev of cluster precursor m/z')
    plt.show()
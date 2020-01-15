from pathlib import Path
import pymzml

import numpy as np


def parse_maracluster_output(path):
    """Parse the output of MaRaCluster.

    Arguments:
        path (str): Path to the data.
    
    Returns:
        A list of numpy arrays, each with indices of clustered spectra.
    """
    path = Path(path)
    clust2indices = []
    clust = []
    with path.open('r') as f:
        for l in f:
            if l == '\n':
                clust2indices.append(clust)
                clust = []
            else:
                clust.append(int(l.split()[1]))
    return [np.array(c) for c in clust2indices]


def get_spectra_array(path2mzml):
    """Get array of spectra.
    
    Arguments:
        path2mzml (str): Path to mzml file.

    Return:
        np.array: an array with spectra.
    """
    return np.array([(s.mz, s.i)
                     for s in pymzml.run.Reader(str(path2mzml))])

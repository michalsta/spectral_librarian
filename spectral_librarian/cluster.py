import numpy as np
import heapq

from spectral_librarian.peaks import PeakCluster, Peak


def peak_distance_clusters(peaks, distance=0.1):
    """Generate peak clusters if successive peaks are closer than 'distance'. 

    Arguments:
        peaks to cluster.

    Yields:
        PeakClusters
    """
    prev_peak = next(peaks)
    peak_cluster = PeakCluster()
    peak_cluster.append(prev_peak)
    for peak in peaks:
        if prev_peak.mz + distance < peak.mz:
            yield peak_cluster
            peak_cluster = PeakCluster()
        peak_cluster.append(peak)
        prev_peak = peak
    yield peak_cluster


def clusters_close(c0, c1, big_mz=1.2):
    """Check, if two clusters are closer than a big jump in m/z."""
    return abs(c1[0].mz - c0[-1].mz) < big_mz 


def filter_far_or_small(clusters, min_el, big_mz=1.2):
    """Filter clusters that are not numerous and separate from other clusters.

    These clusters cannot be part of any isotopic distribution.
    """
    c__ = next(clusters)
    _c_ = next(clusters)
    if clusters_close(c__, _c_, big_mz) or len(c__) >= min_el:
        yield c__
    for __c in clusters:
        if clusters_close(c__, _c_, big_mz) or \
           clusters_close(_c_, __c, big_mz) or \
           len(_c_) >= min_el:
            yield _c_
        c__ = _c_
        _c_ = __c
    if clusters_close(c__, _c_, big_mz) or len(_c_) >= min_el:
        yield _c_



class SpectraCluster(list):
    """Store one MaRaCluster of spectra."""
    def precursor_charges(self):
        return np.array([s.precursor_charge for s in self])

    def precursor_mzs(self):
        return np.array([s.precursor_mz for s in self])

    def precursor_charge(self):
        q = set(self.precursor_charges())
        if len(q) > 1:
            raise RuntimeError("There is no one common precursor charge.")
        return list(q)[0] 

    def __repr__(self):
        g = super().__repr__().replace(', ',',\n')
        return f"SpectraCluster({len(self)})"

    def plot(self, color_map='gist_rainbow', vlines_kwds={}, show=True):
        import matplotlib.pyplot as plt
        cm = plt.get_cmap(color_map)
        for i, spec in enumerate(self):
            col = cm(i/len(self))
            plt.vlines(x=spec.mz, ymin=0, ymax=spec.intensity, colors=col, **vlines_kwds)
        if show:
            plt.show()

    def peaks(self):
        """Iterate over tuples (m/z, intensity, spec_number) for all spectra in a collections.
        Result sorted by m/z."""
        yield from heapq.merge(*([Peak(mz,i,n) for mz,i in zip(s.mz, s.intensity)
                                 ] for n, s in enumerate(self)))

    def refine_intracluster_distance(self, initial_distance=.1, quorum=.9):
        """Refine the initial distance between the peaks.

        To estimate the correct threshold, we first select only peak clusters that are already composed of different peaks.
        Then, we work under the assumption, that close distance between peaks should be more frequent than long jumps.
        If the true number of clusters is K, and number of peaks is N, then N-K-1 first m/z differences in an nondecreasing sequence should correspond to distances between consecutive groups of clusters.
        Here, we plug in the initial estimate of K resulting from the choice of the 'initial_distance'.

        Args:
            initial_distance (nonnegative float): initial distance for peaks.
            quorum (float): Between 0 and 1: fraction of peaks needed to represent the spectrum.
        Returns:
            A refined threshold (minimum from the original and estimated).
        """
        MZ = []
        peak_no = 0
        peak_clusters_no = 0
        for pc in peak_distance_clusters(self.peaks(), initial_distance):
            if len(pc) > len(self)*quorum and pc.peaks_from_different_spectra():
                for p in pc:
                    MZ.append(p.mz)
                peak_no += len(pc)
                peak_clusters_no += 1
        dMZ = np.diff(MZ)
        quantile_distance = np.quantile(dMZ, (peak_no-1-peak_clusters_no)/(peak_no-1))
        if quantile_distance > 0 and quantile_distance < initial_distance/2.0:
            return quantile_distance 
        else:
            return initial_distance/2.0

    def peak_clusters(self, quorum=.5, distance=.1, big_mz=1.2):
        """Generate consecutive clusters of peaks.

        These result from the distance between m/z.
        Small disjoint peaks are not taken into account.
        """
        theshold_clusters = peak_distance_clusters(self.peaks(), distance)
        yield from filter_far_or_small(theshold_clusters, quorum*len(self), big_mz)


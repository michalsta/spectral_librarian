from collections import Counter, OrderedDict, namedtuple


Peak = namedtuple('Peak','mz i spec_no')


class PeakCluster(list):
    def which_spectra(self): 
        return Counter(p.spec_no for p in self)

    def peaks_from_different_spectra(self):
        return all(v==1 for v in self.which_spectra().values())

    def plot(self, vlines_kwds={'color':'orange'}, show=True):
        import matplotlib.pyplot as plt
        MZ, I, _ = list(zip(*self))
        plt.vlines(x=MZ, ymin=0, ymax=I, **vlines_kwds)
        for mz,i,spec_no in self:
            plt.text(x=mz, y=i, s=spec_no)
        if show:
            plt.show()

    def eliminate_less_intense_replicate_peaks(self):
        good_peaks = OrderedDict()
        for mz,i,n in self:
            if n not in good_peaks:
                good_peaks[n] = (mz,i)
            else:
                mz_,i_ = good_peaks[n]
                if i > i_:
                    good_peaks[n] = (mz,i)
        clean_peak_cluster = PeakCluster()
        for n, (mz,i) in good_peaks.items():
            clean_peak_cluster.append((mz,i,n))
        return clean_peak_cluster

    def __repr__(self):
        g = super().__repr__()
        return f"PeakCluster({g[1:-1]})"

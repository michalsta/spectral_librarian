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

    def split(self, mz_split):
        left = PeakCluster()
        right = PeakCluster()
        for p in self:
            if p.mz <= mz_split:
                left.append(p)
            else:
                right.append(p)
        return left, right

    def __repr__(self):
        return f"PeakCluster({len(self)})"

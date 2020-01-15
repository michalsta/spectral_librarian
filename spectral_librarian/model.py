import math


# Py2/3 compat
try:
    range = xrange
except NameError:
    pass



class SpectralModel:
    def __init__(self):
        self._peak_clusters = None

    @staticmethod
    def fromMZPList(L, cluster_gap = 0.01):
        self = SpectralModel()
        self.cluster_gap = cluster_gap

        self.size = len(L)

        spectra = []
        all_confs = []
        for idx, spectrum in enumerate(L):
            masses = spectrum[0]
            probs = spectrum[1]
            spectrum = list(zip(masses, probs))
            norm_factor = math.fsum(x[1] for x in spectrum)
            spectrum = [(x[0], x[1]/norm_factor) for x in spectrum]
            spectra.append(spectrum)

        self.spectra = spectra

        return self

    def get_clusters(self):
        if self._peak_clusters is None:
            self.compute_clusters()
        return self._peak_clusters



    def compute_clusters(self):
        all_confs = []

        for idx, spectrum in enumerate(self.spectra):
            all_confs.extend((x[0], x[1], idx) for x in spectrum)

        all_confs.sort(key = lambda x: x[0])

        peak_clusters = []

        clust_start = 0
        prev_conf = all_confs[0]

        ii = 1
        while ii < len(all_confs):
            curr_conf = all_confs[ii]
            if curr_conf[0] - self.cluster_gap < prev_conf[0]:
                pass
            else:
                peak_clusters.append(all_confs[clust_start:ii])
                clust_start = ii
            prev_conf = curr_conf
            ii += 1

        self._peak_clusters = peak_clusters

    def remove_noise_clusters(self, threshold = 0.2):
        out_clusters = []

        threshold = threshold * self.size

        for cluster in self.get_clusters():
            support = set(x[2] for x in cluster)
            if len(support) > threshold:
                out_clusters.append(cluster)

        self.spectra = [[] for _ in range(self.size)]

        for cluster in out_clusters:
            for mz, prob, sp_id in cluster:
                self.spectra[sp_id].append((mz, prob))

        self.renormalize_spectra()

    def norms(self):
        return [math.fsum(x[1] for x in spectrum) for spectrum in self.spectra]

    def renormalize_spectra(self):
        self.spectra = [[(x[0], x[1]/norm) for x in spectrum] for spectrum, norm in zip(self.spectra, self.norms())]
        self._peak_clusters = None

    def supports(self):
        return [len(set(x[2] for x in spectrum)) for spectrum in self.get_clusters()]

    def plot(self,
             vlines_kwds={},
             span_kwds={'alpha':.1, 'color':'orange'},
             annotate=False,
             show=True):
        import matplotlib.pyplot as plt
        for cl in self.get_clusters():
            min_mz = cl[0][0] - self.cluster_gap
            max_mz = cl[-1][0] + self.cluster_gap
            plt.axvspan(min_mz, max_mz, **span_kwds)
            mz, i, c = zip(*cl)        
            plt.vlines(x=mz, ymin=0, ymax=i, **vlines_kwds)
            if annotate:
                for mz, i, c in cl:
                    plt.text(x=mz, y=i, s=c)
        if show:
            plt.show()



if __name__ == '__main__':
    from collections import Counter
    from pprint import pprint
    import json
    with open("../data/clust0.json") as f:
        S = json.load(f)
    #example = [0.0, 0.0001, 0.0002, 0.0005, 1.0, 1.1, 1.11, 1.111, 1.11111, 1.111111111, 1.112, 1.12, 1.2, 2.0]
    #example = [[(x, 1.0) for x in example]]
    example = S
    SM = SpectralModel.fromMZPList(example)
    C=Counter(SM.supports())
    pprint(C)
    SM.remove_noise_clusters(0.2)
    pprint(SM.get_clusters())
    pprint(SM.supports())
    SM.plot()

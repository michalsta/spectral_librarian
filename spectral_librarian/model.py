import math


# Py2/3 compat
try:
    range = xrange
except NameError:
    pass



class SpectralModel:
    def __init__(self):
        pass
    @staticmethod

    def fromMZPList(L, cluster_gap = 0.01):
        self = SpectralModel()

        spectra = []
        all_confs = []
        for idx, spectrum in enumerate(L):
            print(spectrum)
            norm_factor = math.fsum(x[1] for x in spectrum)
            spectrum = [(x[0], x[1]/norm_factor) for x in spectrum]
            spectra.append(spectrum)
            all_confs.extend((x[0], x[1], idx) for x in spectrum)

        self.spectra = spectra

        all_confs.sort(key = lambda x: x[0])

        peak_clusters = []

        clust_start = 0
        prev_conf = all_confs[0]

        ii = 1
        while ii < len(all_confs):
            curr_conf = all_confs[ii]
            if curr_conf[0] - cluster_gap < prev_conf[0]:
                pass
            else:
                peak_clusters.append(all_confs[clust_start:ii])
                clust_start = ii
            prev_conf = curr_conf
            ii += 1

        self.peak_clusters = peak_clusters

        return self

if __name__ == '__main__':
    example = [0.0, 0.0001, 0.0002, 0.0005, 1.0, 1.1, 1.11, 1.111, 1.11111, 1.111111111, 1.112, 1.12, 1.2, 2.0]
    example = [[(x, 1.0) for x in example]]
    print(SpectralModel.fromMZPList(example).peak_clusters)

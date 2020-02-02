from collections import defaultdict
from math import exp, log, inf
from scipy.spatial import cKDTree as KDTree
from IsoSpecPy import IsoThreshold
from spectral_librarian.cluster import SpectraCluster
from numba import njit, jit, jitclass
import numba.types as nt
import numba

proton_mass = IsoThreshold(0.0, "H1").masses[0]
electron_mass = IsoThreshold(0.0, "E1").masses[0]


from IsoSpecPy.Formulas import *
from IsoSpecPy import *

vect_type = nt.Tuple([nt.double, nt.double, nt.double])
'''@jitclass([('parents', nt.DictType(vect_type, vect_type)),
           ('sizes', nt.DictType(vect_type, nt.int64)),
           ('ex', vect_type)])'''
class DisjointSets(object):
    def __init__(self, v_list : nt.List(vect_type)):
        p = {}
        for x in v_list:
            p[x] = x
        self.parents = p
#        self.parents = {x : x for x in v_list}
        self.sizes = {v_list[0] : 0}
        self.ex = v_list[0]
        del self.sizes[v_list[0]]

    def find(self, what):
        p = self.parents
        while p[what] != what:
            t = p[p[what]]
            p[what] = t
            what = t
        return what

    def union(self, x, y):
        xr = self.find(x)
        yr = self.find(y)
        if xr == yr:
            return
        s = self.sizes
        if s.get(xr, 0) < s.get(yr, 0):
            self.parents[yr] = xr
            s[xr] = s.get(xr, 0) + s.get(yr, 0)
        else:
            self.parents[xr] = yr
            s[yr] = s.get(yr, 0) + s.get(xr, 0)

    def get_sets(self):
        DS = {self.ex : set([self.ex])} #numba.typed.Dict()
        del DS[self.ex]
        f = self.find
        for x in self.parents:
            # Modifying values assigned to keys should be safe while iterating over dict, as long as we don't add/del keys
            fx = f(x)
            if fx in DS:
                DS[fx].append(x)
            else:
                DS[fx] = [x]
        return DS.values()

#@njit 
def find(p, what):
    while p[what] != what:
        t = p[p[what]]
        p[what] = t
        what = t
    return what

#@njit
def disjoint_sets_cluster(pairs):
    parents = {}
    sizes = {}
    for x, y in pairs:
        parents[x] = x
        parents[y] = y
        sizes[x] = 0
        sizes[y] = 0

    for x, y in pairs:
        xr = find(parents, x)
        yr = find(parents, y)

        if xr == yr:
            continue

        s = sizes
        if s.get(xr, 0) < s.get(yr, 0):
            parents[yr] = xr
            s[xr] = s.get(xr, 0) + s.get(yr, 0)
        else:
            parents[xr] = yr
            s[yr] = s.get(yr, 0) + s.get(xr, 0)

    ret = []
    for x in parents:
        ret.append((x, find(parents, x)))
        #ret[x] = find(parents, x)

    return ret
    '''
    ret = {}
    for x in parents:
        rx = find(parents, x)
        if rx in ret:
            ret[rx].append(x)
        else:
            ret[rx] = [x]

    return ret
    '''
        
def get_clusters(vectors, metric, cutoff):
    print("Making KD tree. len(vectors) ==", len(vectors))
    KDT = KDTree(vectors)
    print("KD tree done!")
    pairs = KDT.query_pairs(r=cutoff, p=inf)
    print("pairs done!")
    print("Making DJSet")
#    ds = DisjointSets(vectors)
    print("DJSet done! Making clusters...")
    print("Making list")
    pairs = list(pairs)
    print("calling DSC. len(pairs) ==", len(pairs))
    clusters = disjoint_sets_cluster(pairs)
    print("Clustered, left numba")
    #Actually, gotta invert those. Numba can't for some reason...
    cl = defaultdict(list)
    for pt, idp in clusters:
        cl[idp].append(vectors[pt])
    print("All done")
    return cl.values()
    '''for x1, x2 in pairs:
        ds.union(vectors[x1], vectors[x2])'''
    print("Done!")
    return #ds.get_sets()

    
    '''
    while len(unvisited) > 0:
        print(len(unvisited))
        current_cluster = set()
        dfs_grey = set([unvisited.pop()])
        while len(dfs_grey) > 0:
            current = dfs_grey.pop()
            print(current)
            unvisited.discard(current)
            current_cluster.add(current)
            for x in set():
                if metric(x, current) < cutoff:
                    print("BLA")
                    unvisited.discard(x)
                    dfs_grey.add(x)
        yield current_cluster
'''

#import fastcluster

#def get_clusters(vectors, metric, cutoff):
#    fastcluster.linkage_vector(vectors, metric='chebychev')

class TestCase(SpectraCluster):
    def __init__(self):
        import json
        from pkg_resources import resource_stream
        with resource_stream("spectral_librarian", "data/clust0.json") as f:
            spectra = json.load(f)
        super().__init__(spectra)


def metric_transformation(m1, m2, int1, int2):
#    print(m1, m2, int1, int2)
    ret = (m1, m2, log(int1/int2)/1000.0)
#    print(ret)
    return ret


def sup_metric(e1, e2):
    return max(abs(c1-c2) for c1, c2 in zip(e1, e2))

class EdgesGraph(object):
    def __init__(self, spectra_cluster, charge_mass = proton_mass):
        self.spectra = spectra_cluster
        self.charge_mass = charge_mass

        self.edge_to_spectrum = defaultdict(list)
        edges = []

        for idx, (masses, intensities) in enumerate(self.spectra):
            print(len(masses))
            for m1, int1 in zip(masses, intensities):
                for m2, int2 in zip(masses, intensities):
                    if m2 <= m1: continue
                
                    edge = metric_transformation(m1, m2, int1, int2)
                    self.edge_to_spectrum[edge].append(idx)

                    edges.append(edge)

            print(list(self.edge_to_spectrum.values()).count([idx]))
#            break # for debugging only

        for cluster in get_clusters(edges, sup_metric, 1.0):
            pass
#            print(cluster)
        #kd_tree = KDTree(kd_tree_tbl)
        print("Construction done!")




        

if __name__ == "__main__":
    EdgesGraph(TestCase())

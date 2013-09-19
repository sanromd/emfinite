import cPickle as pickle
import numpy as np
from mpi4py import MPI

class Probe(object):

    def __init__(self, da, entries):
        entries = np.asarray(entries, dtype='i')
        assert entries.ndim == 2
        assert da.dim == entries.shape[1]
        #
        local = []
        (xi, xf), (yi, yf) = da.getRanges()
        for i, j in entries:
            if (xi <= i < xf and
                yi <= j < yf):
                local.append([i, j])
        #
        self.comm = da.comm
        self.da = da
        self.indices = np.asarray(local, dtype='i')
        self.cache = {}

    def probe(self, key, Q):
        entry = self.cache.setdefault(key, {})
        q = self.da.getVecArray(Q)
        for i, j in self.indices:
            i, j = int(i), int(j)
            data = entry.setdefault((i,j), [])
            data.append(q[i,j])

    def save(self, filename):
        comm = MPI.COMM_WORLD
        caches = comm.gather(self.cache, root=0)
        if comm.rank == 0:
            #
            result = caches[0]
            for i in range(1,len(caches)):
                entry = caches[i]
                for key in entry:
                    result[key].update(entry[key])
            #
            fh = open(filename, "wb")
            pickle.dump(result, fh, -1)
            fh.close()


# import cPickle as pickle
# fh = open("probe.dat", "r")
# data = pickle.load(fh)
# fh.close()

import numpy as np

class Reader(object):

    VEC_ID = 1211214
    precision = {
        'single' : {'real' : '>f4', 'complex' : '>c8' },
        'double' : {'real' : '>f8', 'complex' : '>c16'},
        }
    indices = {
        '32bit' : '>i4',
        '64bit' : '>i8',
        }

    def __init__(self,
                 precision='double',
                 scalar='real',
                 indices='32bit'):
        I = self.indices[indices]
        R = self.precision[precision]['real']
        S = self.precision[precision][scalar]
        self._types = tuple(np.dtype(t) for t in (I, R, S))

    def _read(self, fid, dtype, count=None):
        if count is None:
            array = np.fromfile(fid, dtype, 1)[0]
        else:
            array = np.fromfile(fid, dtype, count)
        return array.astype(dtype.newbyteorder('='))

    def _read_vec(self, fh):
        VEC_ID = self.VEC_ID
        I,R,S  = self._types
        _read  = self._read
        #
        clsid = _read(fh, I)
        assert clsid == VEC_ID
        n = _read(fh, I)
        A = _read(fh, S, n)
        assert len(A) == n
        #
        return A
        
    def read(self, filename):
        fh = open(filename, 'rb')
        try:
            Q1 = self._read_vec(fh)
            Q2 = self._read_vec(fh)
            Q3 = self._read_vec(fh)
        finally:
            fh.close()
        return Q1, Q2, Q3


if __name__ == "__main__":
    import sys
    path = sys.argv[1]
    reader = Reader()
    Q1, Q2, Q3 = reader.read(filename)

    from matplotlib import pylab
    nx, ny = 35, 35
    x, y = np.mgrid[0:1:1j*nx, 0:1:1j*ny]
    q1 = Q1.reshape([nx,ny], order='F')
    q2 = Q2.reshape([nx,ny], order='F')
    q3 = Q3.reshape([nx,ny], order='F')
    
    pylab.figure()
    pylab.contourf(x,y,q1)
    pylab.figure()
    pylab.contourf(x,y,q2)
    pylab.figure()
    pylab.contourf(x,y,q3)
    pylab.show()

import numpy as np

from common import HasInner

class NoPrec(HasInner):
    def __init__(self, *args):
        pass

    def solve(self, r):
        return r[self.inner]

class JacDiag(NoPrec):
    def __init__(self, Di):

        #self.Mi = ne.evaluate("""
        self.Mi = 1.0/Di[self.inner]
        #""")

    def solve(self, r):
        #Mi=self.Mi
        #ri = r[self.inner]
        #return ne.evaluate("Mi*r")
        z = self.Mi*r[self.inner]
        return z

class JacBlock(NoPrec):
    def __init__(self, Kx, Ky):
        Kx_c = Kx[self.inner]
        Kx_r = Kx[self.xr_idx]

        Ky_c = Ky[self.inner]
        Ky_r = Ky[self.yr_idx]

        self.a = -Ky_c.copy()
        self.b = (1.0 + (Kx_r + Kx_c) + (Ky_r + Ky_c)).copy()
        self.c = -Ky_r.copy()

        self.stride = 4
        srange = [i for i in xrange(self.stride)]

        self.arrshape = (self.a.size/self.stride, self.stride)

        self.a = self.a.reshape(self.arrshape).copy(order='F')
        self.b = self.b.reshape(self.arrshape).copy(order='F')
        self.c = self.c.reshape(self.arrshape).copy(order='F')

        self.cprime = np.zeros_like(self.c)
        self.bfprime = np.zeros_like(self.c)
        sidx = np.s_[:,0]

        self.cprime[sidx] = self.c[sidx]/self.b[sidx]

        for idx in srange[1:]:
            iidx = np.s_[:,idx]
            midx = np.s_[:,idx-1]

            self.cprime[iidx] = self.c[iidx]/(self.b[iidx] - self.a[iidx]*self.cprime[midx])

            self.bfprime[iidx] = self.b[iidx] - self.a[iidx]*self.cprime[midx]

    def solve(self, r):
        d = r[self.inner]

        stride = self.stride

        d = d.reshape(self.arrshape).copy(order='F')
        x = np.zeros_like(self.c)

        dprime = np.zeros_like(d)

        sidx = np.s_[:,0]
        dprime[sidx] = d[sidx]/self.b[sidx]

        srange = [i for i in xrange(self.stride)]
        for idx in srange[1:]:
            iidx = np.s_[:,idx]
            midx = np.s_[:,idx-1]

            dprime[iidx] = (d[iidx] - self.a[iidx]*dprime[midx]) / self.bfprime[iidx]

        midx = np.s_[:,stride-1]
        x[midx] = dprime[midx]

        for idx in srange[:-1][::-1]:
            iidx = np.s_[:,idx]
            pidx = np.s_[:,idx+1]
            x[iidx] = (dprime[iidx] - self.cprime[iidx]*x[pidx])

        x = x.reshape(r[self.inner].shape)

        return x

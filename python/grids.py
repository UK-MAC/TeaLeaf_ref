import numpy as np
import functools
import scipy as sp
from scipy.linalg.blas import ddot

import numexpr as ne

from preconditioners import *

from common import HasInner, calc_eigs, calc_ch_coefs, pltim

class Grid(HasInner):
    def __init__(self, params, dims):
        self.params = params

        self.Ky = np.zeros((self.params.x_cells+2, self.params.y_cells+2))
        self.Kx = self.Ky.copy()
        self.Di = self.Ky.copy()
        self.density = np.zeros_like(self.Kx)
        self.energy = np.zeros_like(self.Kx)

        dens_in = self.density[self.inner]
        ener_in = self.energy[self.inner]

        dx = ((self.params.xmax - self.params.xmin)/self.params.x_cells)
        dy = ((self.params.ymax - self.params.ymin)/self.params.y_cells)
        self.rx = self.params.dtinit/(dx**2)
        self.ry = self.params.dtinit/(dy**2)

        for state in self.params.states:
            state_info = self.params.states[state]

            if state_info['geometry'] == 0:
                dens_in[:,:] = state_info['density']
                ener_in[:,:] = state_info['energy']
                continue

            xmin = np.rint(state_info['xmin']*self.params.x_cells/self.params.xmax)
            ymin = np.rint(state_info['ymin']*self.params.y_cells/self.params.ymax)
            xmax = np.rint(state_info['xmax']*self.params.x_cells/self.params.xmax)
            ymax = np.rint(state_info['ymax']*self.params.y_cells/self.params.ymax)

            dens_in[xmin:xmax, ymin:ymax] = state_info['density']
            ener_in[xmin:xmax, ymin:ymax] = state_info['energy']

        self.zero_boundaries(self.density)
        self.zero_boundaries(self.energy)

        self.u0 = self.energy*self.density
        self.u = self.u0.copy()

        self.makeA(self.density)

        if self.params.tl_use_dpcg:
            self.dims = dims

            self.x_tiles = self.dims[0]
            self.y_tiles = self.dims[1]

            self.Ky_small = np.zeros((self.dims[0]+2, self.dims[1]+2))
            self.Kx_small = self.Ky_small.copy()
            self.Di_small = self.Ky_small.copy()

            # used kind of like meshgrid()
            self.tile_dims_x = np.zeros(self.dims[0])
            self.tile_dims_y = np.zeros(self.dims[1])

            sub_arrays = np.array_split(dens_in, self.dims[0], 0)
            for i,s in enumerate(sub_arrays):
                self.tile_dims_x[i] = s.shape[0]

            sub_arrays = np.array_split(dens_in, self.dims[1], 1)
            for i,s in enumerate(sub_arrays):
                self.tile_dims_y[i] = s.shape[1]

            self.x_cumsum = np.cumsum(self.tile_dims_x).astype(np.int32)
            self.x_cumsum_minus = (self.x_cumsum - self.tile_dims_x).astype(np.int32)
            self.y_cumsum = np.cumsum(self.tile_dims_y).astype(np.int32)
            self.y_cumsum_minus = (self.y_cumsum - self.tile_dims_y).astype(np.int32)

            self.makeA_small()

            self.calculated_coarse_eigs = False

        if self.params.tl_preconditioner_type == "none":
            self.M = NoPrec(self.Kx, self.Ky)
        elif self.params.tl_preconditioner_type == "jac_diag":
            self.M = JacDiag(self.Di)
        elif self.params.tl_preconditioner_type == "jac_block":
            self.M = JacBlock(self.Kx, self.Ky)

        self.z, self.p, self.r, self.w = self.init_arrays(self.u)

    def makeA(self, d):
        self.Kx[:] = 0
        self.Ky[:] = 0
        self.Di[:] = 0

        self.Kx[self.xr_idx] = (d[self.inner] + d[self.xr_idx])/(2.0*d[self.inner]*d[self.xr_idx])
        self.Ky[self.yr_idx] = (d[self.inner] + d[self.yr_idx])/(2.0*d[self.inner]*d[self.yr_idx])

        self.zero_boundaries(self.Kx)
        self.zero_boundaries(self.Ky)

        self.Kx *= self.rx
        self.Ky *= self.ry

        Kx_r = self.Kx[self.xr_idx]
        Kx_c = self.Kx[self.inner]

        Ky_r = self.Ky[self.yr_idx]
        Ky_c = self.Ky[self.inner]

        down  = Ky_c
        up    = Ky_r
        left  = Kx_c
        right = Kx_r

        self.Di[self.inner] = 1.0 + (left + right) + (down + up)

    def makeA_small(self):
        self.Kx_small[:] = 0
        self.Ky_small[:] = 0
        self.Di_small[:] = 0

        Kx_inner = self.Kx[self.inner]
        Ky_inner = self.Ky[self.inner]

        Kx_small_inner = self.Kx_small[self.inner]
        Ky_small_inner = self.Ky_small[self.inner]
        Di_small_inner = self.Di_small[self.inner]

        # TODO make numpy
        for i in xrange(self.dims[0]):
            for j in xrange(self.dims[1]):
                Kx_c =   np.sum(Kx_inner[self.x_cumsum_minus[i  ]:self.x_cumsum[i  ], self.y_cumsum_minus[j  ]:self.y_cumsum[j  ]])

                try:
                    Kx_r=np.sum(Kx_inner[self.x_cumsum_minus[i+1]:self.x_cumsum[i+1], self.y_cumsum_minus[j  ]:self.y_cumsum[j  ]])
                except IndexError as e:
                    Kx_r=0

                Ky_c =   np.sum(Ky_inner[self.x_cumsum_minus[i  ]:self.x_cumsum[i  ], self.y_cumsum_minus[j  ]:self.y_cumsum[j  ]])

                try:
                    Ky_r=np.sum(Ky_inner[self.x_cumsum_minus[i  ]:self.x_cumsum[i  ], self.y_cumsum_minus[j+1]:self.y_cumsum[j+1]])
                except IndexError as e:
                    Ky_r=0

                if i == 0:
                    Kx_c = 0
                if j == 0:
                    Ky_c = 0

                Kx_small_inner[i  ,j  ]=Kx_c
                Ky_small_inner[i,  j  ]=Ky_c

                #if i < self.dims[0] - 1:
                #    Kx_small_inner[i+1,j  ]=Kx_r
                #if j < self.dims[1] - 1:
                #    Ky_small_inner[i,  j+1]=Ky_r

                Di_small_inner[i,j]  = (self.x_cumsum[i]-self.x_cumsum_minus[i])*(self.y_cumsum[j]-self.y_cumsum_minus[j]) # sum the diagonals
                Di_small_inner[i,j] += (Kx_c + Kx_r) + (Ky_c + Ky_r)

        self.zero_boundaries(self.Kx_small)
        self.zero_boundaries(self.Ky_small)

    def dotvec(self, a, b):
        an = a[self.inner].ravel()
        bn = b[self.inner].ravel()

        #norm = ne.evaluate("""sum(an*bn)""")

        norm = ddot(an, bn)

        return norm

    def init_arrays(self, u):
        p = np.zeros_like(u)
        r = np.zeros_like(u)
        w = np.zeros_like(u)
        z = np.zeros_like(u)

        w[self.inner] = self.matmul(u)
        r[self.inner] = self.u0[self.inner] - w[self.inner]

        self.zero_boundaries(z)
        self.zero_boundaries(p)
        self.zero_boundaries(r)
        self.zero_boundaries(w)

        return z, p, r, w

    def calc_residual(self, w, r, b, x, matmul_func):
        # r = b - Ax
        w[self.inner] = matmul_func(x)
        r[self.inner] = b[self.inner] - w[self.inner]

    def exactrro(self, b, x, matmul_func):
        r = np.zeros_like(b)
        w = np.zeros_like(b)
        self.calc_residual(w, r, b, x, matmul_func)
        return self.dotvec(r, r)

    def zero_boundaries(self, arr):
        outer_cols_idx = np.s_[:,(0,-1)]
        outer_rows_idx = np.s_[(0,-1),:]

        arr[outer_rows_idx] = 0
        arr[outer_cols_idx] = 0

    def update_halo(self, *args):
        if self.params.reflective_boundary:
            for array in args:
                array[ 0,:] = array[ 1,:]
                array[-1,:] = array[-2,:]
                array[:, 0] = array[:, 1]
                array[:,-1] = array[:,-2]
        else:
            # possible overkill
            for array in args:
                self.zero_boundaries(array)

    def ppcg(self, u, p, r, sd, w, z, theta, alphas, betas, inner_its, rro, matmul_func, M=NoPrec(0,0)):
        w[self.inner] = matmul_func(p)
        pw = self.dotvec(p, w)

        alpha = rro/pw

        u += alpha*p
        r -= alpha*w

        z[self.inner] = M.solve(r)
        sd = z*(1.0/theta)

        for i in xrange(inner_its):
            w[self.inner] = matmul_func(sd)
            r -= w
            u += sd

            z[self.inner] = M.solve(r)

            #sd = alphas[i]*sd + betas[i]*z
            sd *= alphas[i]
            sd += betas[i]*z

        rrn = self.dotvec(r, z)

        beta = rrn/rro

        #p[:] = z + beta*p
        p *= beta
        p += z

        return rrn

    def cg(self, u, p, r, w, z, rro, matmul_func, M=NoPrec(0,0)):
        w[self.inner] = matmul_func(p)
        pw = self.dotvec(p, w)

        alpha = rro/pw

        u += alpha*p
        r -= alpha*w

        z[self.inner] = M.solve(r)
        rrn = self.dotvec(r, z)

        beta = rrn/rro

        #p[:] = z + beta*p
        p *= beta
        p += z

        self.update_halo(u, p)

        return rrn, alpha, beta

    def smooth(self, arr):
        ux_r = arr[self.xr_idx]
        ux_l = arr[self.xl_idx]

        uy_r = arr[self.yr_idx]
        uy_l = arr[self.yl_idx]

        ux_c = arr[self.inner]

        Kx_r = self.Kx[self.xr_idx]
        Kx_c = self.Kx[self.inner]

        Ky_r = self.Ky[self.yr_idx]
        Ky_c = self.Ky[self.inner]

        diag = self.Di[self.inner]

        ins = arr[self.inner]

        un = ins.copy()

        u0 = self.u0[self.inner]

        damp = 2.0/3.0

        #arr[self.inner] = ne.evaluate("""damp*(u0  + ((Kx_r*ux_r + Kx_c*ux_l) + (Ky_r*uy_r + Ky_c*uy_l)))\
        #        /diag + (1 - damp)*ux_c""")

        arr[self.inner] = damp*(u0  + ((Kx_r*ux_r + Kx_c*ux_l) + (Ky_r*uy_r + Ky_c*uy_l)))\
                /diag + (1 - damp)*ux_c

    def matmul(self, arr):
        ux_r = arr[self.xr_idx]
        ux_l = arr[self.xl_idx]

        ins = arr[self.inner]

        uy_r = arr[self.yr_idx]
        uy_l = arr[self.yl_idx]

        Kx_r = self.Kx[self.xr_idx]
        Kx_c = self.Kx[self.inner]

        Ky_r = self.Ky[self.yr_idx]
        Ky_c = self.Ky[self.inner]

        diag = self.Di[self.inner]

        #return ne.evaluate("""diag*ins - (Kx_r*ux_r + Kx_c*ux_l) - (Ky_r*uy_r + Ky_c*uy_l)""")
        return diag*ins - (Kx_r*ux_r + Kx_c*ux_l) - (Ky_r*uy_r + Ky_c*uy_l)

    def matmul_small(self, arr):
        ux_r = arr[self.xr_idx]
        ux_l = arr[self.xl_idx]

        ins = arr[self.inner]

        uy_r = arr[self.yr_idx]
        uy_l = arr[self.yl_idx]

        Kx_r = self.Kx_small[self.xr_idx]
        Kx_c = self.Kx_small[self.inner]

        Ky_r = self.Ky_small[self.yr_idx]
        Ky_c = self.Ky_small[self.inner]

        diag = self.Di_small[self.inner]

        #return ne.evaluate("""diag*ins - (Kx_r*ux_r + Kx_c*ux_l) - (Ky_r*uy_r + Ky_c*uy_l)""")
        return diag*ins - (Kx_r*ux_r + Kx_c*ux_l) - (Ky_r*uy_r + Ky_c*uy_l)

    def iter_sub_array(self, big_array, small_array, function):
        def get_slice_and_call(((x,y),_)):
            arr_slice = np.s_[self.x_cumsum_minus[x]:self.x_cumsum[x], self.y_cumsum_minus[y]:self.y_cumsum[y]]
            function(x, y, big_array, small_array, arr_slice)

        map(get_slice_and_call, np.ndenumerate(small_array[self.inner]))

    def restrict_Zt(self, r):
        t2 = np.zeros((2+self.dims[0], 2+self.dims[1]))

        def sum_sub(x, y, big_array, small_array, arr_slice):
            small_array[self.inner][x, y] = np.sum(big_array[self.inner][arr_slice])
            #if not x and not y: print small_array[self.inner][x, y] , arr_slice

        self.iter_sub_array(r, t2, sum_sub)

        return t2

    def prolong_Z(self, t2):
        prolonged = np.zeros_like(self.u[self.inner])

        def broadcast(x, y, big_array, small_array, arr_slice):
            big_array[arr_slice][:] = small_array[self.inner][x, y]

        self.iter_sub_array(prolonged, t2, broadcast)

        return prolonged

    def matmul_ZtA(self, z):
        t1 = np.zeros((2+self.dims[0], 2+self.dims[1]))

        def sum_and_multiply(x, y, big_array, small_array, arr_slice):
            small_array[self.inner][x, y] = np.sum(big_array[arr_slice])

        self.iter_sub_array(self.matmul(z), t1, sum_and_multiply)

        return t1

    def solve_E(self, x_s, b_s):
        r_s = np.zeros_like(x_s)
        w_s = np.zeros_like(x_s)
        p_s = np.zeros_like(x_s)
        z_s = np.zeros_like(x_s)
        sd_s = np.zeros_like(x_s)

        matmul_bound = functools.partial(self.matmul_small)

        #M_s = NoPrec(E, E)
        #M_s = JacDiag(self.Di_small)
        if self.params.tl_preconditioner_type == "jac_diag":
            M_s = JacDiag(self.Di_small)
        elif self.params.tl_preconditioner_type == "jac_block":
            M_s = JacBlock(self.Kx_small, self.Ky_small)
        else:
            M_s = NoPrec(self.Kx, self.Ky)

        # use zero initial guess which allows for inexact coarse solves
        x_s[:] = 0

        self.calc_residual(w_s, r_s, b_s, x_s, matmul_func=matmul_bound)

        z_s[self.inner] = M_s.solve(r_s)
        p_s[:] = z_s

        rro = self.dotvec(r_s, p_s)

        initial = rro

        inner_its = 4

        cg_alphas = []
        cg_betas = []

        i = 0

        for i in xrange(self.params.coarse_solve_max_iters):
            if 1:
                rro, alpha, beta = self.cg(x_s, p_s, r_s, w_s, z_s, rro, self.matmul_small, M=M_s)
            else:
                if self.calculated_coarse_eigs:
                    rro = self.ppcg(x_s, p_s, r_s, sd_s, w_s, z_s,
                        self.inner_theta, self.inner_ch_alphas, self.inner_ch_betas,
                        inner_its, rro, matmul_bound, M=M_s)
                else:
                    eigmin = 0.01
                    eigmax = 2

                    # use above bounds if diagonal jacobi is used
                    if 1:
                        for j in xrange(35):
                            rro, alpha, beta = self.cg(x_s, p_s, r_s, w_s, z_s, rro, matmul_bound, M=M_s)
                            if rro == 0:
                                break

                            cg_alphas.append(alpha)
                            cg_betas.append(beta)

                        i += j

                        eigmin, eigmax = calc_eigs(cg_alphas, cg_betas)

                    self.inner_theta, self.inner_ch_alphas, self.inner_ch_betas = calc_ch_coefs(eigmin, eigmax,
                        self.dotvec(self.u0, self.u0),
                        self.dotvec(self.r, self.r),
                        self.params.coarse_solve_eps
                        )

                    #self.calculated_coarse_eigs = True

                    #print eigmin, eigmax, eigmax/eigmin

            if np.sqrt(abs(rro)) < self.params.coarse_solve_eps*np.sqrt(abs(initial)):
            #if np.sqrt(abs(rro)) < 0.9*np.sqrt(abs(initial)):
                break

        #print initial
        #print i, "iterations", rro

        #self.calc_residual(w_s, r_s, b_s, x_s, matmul_func=matmul_bound)
        #print initial, i, rro, self.exactrro(b_s, x_s, matmul_bound)

        return x_s

    def solve_coarse(self, r, z):
        # 12/13
        z[self.inner] = self.M.solve(r)

        # 14
        t1 = self.matmul_ZtA(z)

        # 15
        t2 = self.restrict_Zt(r)

        # 16
        t1 -= t2

        # 17
        t2 = self.solve_E(t2, t1)

        # 18/19
        z[self.inner] -= self.prolong_Z(t2)

    def init_fadef(self, u, p, r, w, z):
        # 4/5
        self.calc_residual(w, r, self.u0, u, matmul_func=self.matmul)

        if self.params.tl_use_dpcg:
            # 6
            t1 = self.restrict_Zt(r)

            # 7
            t2 = -t1

            t2 = self.solve_E(t2, t1)

            # 8/9
            u[self.inner] -= self.prolong_Z(t2)

            # 10/11
            self.calc_residual(w, r, self.u0, u, matmul_func=self.matmul)

            # 12-19
            self.solve_coarse(r, z)

        elif self.params.tl_use_cg:
            z[self.inner] = self.M.solve(r)

        # 20
        p[:] = z

    def fadef(self, u, p, r, w, z):
        # 2
        w[self.inner] = self.matmul(p)

        # 3
        r_old = r.copy()
        rro = self.dotvec(r, z)
        pw = self.dotvec(p, w)
        alpha = rro/pw

        # 4
        u += alpha*p

        # 5
        r -= alpha*w

        # 6-13
        self.solve_coarse(r, z)

        # 14
        rrn = self.dotvec(r - r_old, z)
        beta = rrn/rro

        # 15
        p[self.inner] = z[self.inner] + beta*p[self.inner]

        return rrn, alpha, beta


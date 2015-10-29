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
        self.density = np.zeros_like(self.Kx)
        self.energy = np.zeros_like(self.Kx)

        dx = ((self.params.xmax - self.params.xmin)/self.params.x_cells)
        dy = ((self.params.ymax - self.params.ymin)/self.params.y_cells)

        self.rx = self.params.dtinit/(dx**2)
        self.ry = self.params.dtinit/(dy**2)

        self.dims = dims

        self.x_tiles = self.dims[0]
        self.y_tiles = self.dims[1]

        dens_in = self.density[self.inner]
        ener_in = self.energy[self.inner]

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

        if self.params.tl_preconditioner_type == "none":
            self.M = NoPrec(self.Kx, self.Ky)
        elif self.params.tl_preconditioner_type == "jac_diag":
            self.M = JacDiag(self.Kx, self.Ky)
        elif self.params.tl_preconditioner_type == "jac_block":
            self.M = JacBlock(self.Kx, self.Ky)

        self.z, self.p, self.r, self.w = self.init_arrays(self.u)

        self.calculated_coarse_eigs = False

    def makeA(self, d):
        self.Kx[:] = 0
        self.Ky[:] = 0

        self.Kx[self.xr_idx] = (d[self.inner] + d[self.xr_idx])/(2.0*d[self.inner]*d[self.xr_idx])
        self.Ky[self.yr_idx] = (d[self.inner] + d[self.yr_idx])/(2.0*d[self.inner]*d[self.yr_idx])

        self.zero_boundaries(self.Kx)
        self.zero_boundaries(self.Ky)

        self.Kx *= self.rx
        self.Ky *= self.ry

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
        r[:] = b - w

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

    def sum_matrix_row(self):
        Kx_r = self.Kx[self.xr_idx]
        Kx_c = self.Kx[self.inner]

        Ky_r = self.Ky[self.yr_idx]
        Ky_c = self.Ky[self.inner]

        down = -Ky_c
        up = -Ky_r
        left = -Kx_c
        right = -Kx_r

        sums = ne.evaluate("(1.0 - down - up - left - right) + down + up + left + right")

        return sums

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

        ins = arr[self.inner]

        un = ins.copy()

        u0 = self.u0[self.inner]

        damp = 2.0/3.0

        arr[self.inner] = ne.evaluate("""damp*(u0  + ((Kx_r*ux_r + Kx_c*ux_l) + (Ky_r*uy_r + Ky_c*uy_l)))\
                /(1.0 + (Kx_c+Kx_r) + (Ky_c+Ky_r))\
                + (1 - damp)*ux_c""")

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

        #return ne.evaluate("""
        return (1.0 + (Kx_r + Kx_c)
                + (Ky_r + Ky_c))*ins             \
                - (Kx_r*ux_r + Kx_c*ux_l)        \
                - (Ky_r*uy_r + Ky_c*uy_l)
        #""")

    def iter_sub_array(self, big_array, small_array, function):
        def get_slice_and_call(((x,y),_)):
            arr_slice = np.s_[self.x_cumsum_minus[x]:self.x_cumsum[x], self.y_cumsum_minus[y]:self.y_cumsum[y]]
            function(x, y, big_array, small_array, arr_slice)

        #for (x,y),_ in np.ndenumerate(small_array[self.inner]):
        #    get_slice_and_call((x,y),_)

        map(get_slice_and_call, np.ndenumerate(small_array[self.inner]))

    def restrict_Zt(self, arr):
        t2 = np.zeros((2+self.dims[0], 2+self.dims[1]))

        def sum_sub(x, y, big_array, small_array, arr_slice):
            small_array[self.inner][x, y] = np.sum(big_array[self.inner][arr_slice])
            #if not x and not y: print small_array[self.inner][x, y] , arr_slice

        self.iter_sub_array(arr, t2, sum_sub)

        return t2

    def prolong_Z(self, small_arr):
        prolonged = np.zeros_like(self.u[self.inner])

        def broadcast(x, y, big_array, small_array, arr_slice):
            big_array[self.inner][arr_slice][:] = small_array[self.inner][x, y]

        self.iter_sub_array(prolonged, small_arr, broadcast)

        return prolonged

    def matmul_ZtA(self, arr):
        t1 = np.zeros((2+self.dims[0], 2+self.dims[1]))

        def sum_and_multiply(x, y, big_array, small_array, arr_slice):
            small_array[self.inner][x, y] = np.sum(big_array[self.inner][arr_slice])

        self.iter_sub_array(self.matrix_row_sums*arr[self.inner], t1, sum_and_multiply)

        return t1

    def solve_E(self, x_s, b_s):
        r_s = np.zeros_like(x_s)
        w_s = np.zeros_like(x_s)
        p_s = np.zeros_like(x_s)
        z_s = np.zeros_like(x_s)
        sd_s = np.zeros_like(x_s)

        E = np.zeros_like(x_s)

        def make_e(x, y, big_array, small_array, arr_slice):
            small_array[self.inner][x, y] = np.sum(big_array[arr_slice])

        self.iter_sub_array(self.matrix_row_sums, E, make_e)

        def matmul_small(arr, E):
            ux_r = arr[self.xr_idx]
            ux_l = arr[self.xl_idx]

            ins = arr[self.inner]

            uy_r = arr[self.yr_idx]
            uy_l = arr[self.yl_idx]

            Kx_r = E[self.xr_idx]
            Kx_c = E[self.inner]

            Ky_r = E[self.yr_idx]
            Ky_c = E[self.inner]

            #return ne.evaluate("""
            return (1.0 + (Kx_r + Kx_c)
                    + (Ky_r + Ky_c))*ins             \
                    - (Kx_r*ux_r + Kx_c*ux_l)        \
                    - (Ky_r*uy_r + Ky_c*uy_l)
            #""")

        matmul_bound = functools.partial(matmul_small, E=E)

        #M_s = NoPrec(E, E)
        M_s = JacDiag(E, E)

        self.calc_residual(w_s, r_s, b_s, x_s, matmul_func=matmul_bound)

        z_s[self.inner] = M_s.solve(r_s)
        p_s[:] = z_s

        rro = self.dotvec(r_s, p_s)

        cg_alphas = []
        cg_betas = []

        i = 0

        initial = rro

        inner_its = 4

        for i in xrange(100):
            if self.calculated_coarse_eigs:
                rro = self.ppcg(x_s, p_s, r_s, sd_s, w_s, z_s,
                    self.inner_theta, self.inner_ch_alphas, self.inner_ch_betas,
                    inner_its, rro, matmul_bound, M=M_s)
            else:
                eigmin = 0.01
                eigmax = 2

                # use above bounds if diagonal jacobi is used
                if 0:
                    for i in xrange(35):
                        rro, alpha, beta = self.cg(x_s, p_s, r_s, w_s, z_s, rro, matmul_bound, M=M_s)
                        if rro == 0:
                            break

                        cg_alphas.append(alpha)
                        cg_betas.append(beta)

                    eigmin, eigmax = calc_eigs(cg_alphas, cg_betas)

                self.inner_theta, self.inner_ch_alphas, self.inner_ch_betas = calc_ch_coefs(eigmin, eigmax, 
                    self.dotvec(self.u0, self.u0),
                    self.dotvec(self.r, self.r),
                    #np.finfo(type(self.r[0,0])),
                    self.params.tl_eps
                    )
                self.calculated_coarse_eigs = True

            if np.sqrt(abs(rro)) < self.params.tl_eps*np.sqrt(abs(initial)):
                break

        return x_s

    def solve_coarse(self, r, z):
        # 12/13
        z[self.inner] = self.M.solve(r)

        # 14
        t1 = self.matmul_ZtA(z)

        # 15
        t2 = self.restrict_Zt(r)

        # 16
        t1 += t2

        # 17
        t2 = self.solve_E(t2, t1)

        # 18/19
        z[self.inner] -= self.prolong_Z(t2)

    def init_fadef(self, u, p, r, w, z):
        # 4/5
        self.calc_residual(w, r, self.u0, u, matmul_func=self.matmul)

        if self.params.tl_use_dpcg:
            # -1
            self.matrix_row_sums = self.sum_matrix_row()

            # 6
            t2 = self.restrict_Zt(r)

            # 7
            t1 = t2.copy()
            t2 = self.solve_E(t2, t1)

            # 8/9
            u[self.inner] += self.prolong_Z(t2)

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
        # TODO find out why it doesn't converge without this
        z[self.inner] = self.M.solve(r)
        p[self.inner] = z[self.inner] + beta*p[self.inner]

        return rrn, alpha, beta


#!/usr/bin/env python

# TODO
# parallel with multiprocessing

import itertools
import numpy as np

from grids import Grid

from ctypes import cdll, c_int, byref

import time

from tea_params import Params

from common import calc_eigs, calc_ch_coefs

np.set_printoptions(linewidth=1000, precision=4)

class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

class Solver(object):
    def __init__(self):
        self.params = Params("../tea.in")

        libmpi = cdll.LoadLibrary("libmpi.so")
        libmpi.MPI_Init(None, None)
        nnodes = c_int(self.params.tiles_per_task)
        ndims = c_int(2)
        dims = (c_int*2)()
        dims[0] = 0
        dims[1] = 0

        libmpi.MPI_Dims_create(nnodes, ndims, byref(dims))

        dims = dims[::-1]

        libmpi.MPI_Finalize()

        self.figidx = 1

        self.grid = Grid(self.params, dims)

    def animate(self):
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
        from matplotlib.colors import LogNorm

        def getarr():
            return self.grid.u
        fig, ax = plt.subplots(1, 1)

        im = ax.imshow(getarr(),
            cmap=plt.get_cmap('seismic'),
            #norm=LogNorm(vmin=0.1, vmax=np.max(getarr())),
            )

        for (x,y) in itertools.product(enumerate(np.cumsum(self.grid.tile_dims_x)), enumerate(np.cumsum(self.grid.tile_dims_y))):
            print x, y
            ax.axhline(1)

        for x in np.cumsum(self.grid.tile_dims_x):
            ax.axhline(x, c='y')
        for y in np.cumsum(self.grid.tile_dims_y):
            ax.axvline(y, c='y')

        def updatefig(*args):
            i = args[0]

            if self.params.tl_use_cg:
                rro = self.grid.dotvec(self.grid.r, self.grid.z)
                rro, alpha, beta = self.grid.cg(self.grid.u, self.grid.p, self.grid.r, self.grid.w, self.grid.z,
                    rro, self.grid.matmul, M=self.grid.M)
            elif self.params.tl_use_dpcg:
                rro, alpha, beta = self.grid.fadef(self.grid.u, self.grid.p, self.grid.r, self.grid.w, self.grid.z)

            error = np.sqrt(abs(rro))

            if not i % 10:
                print "{0: 4d} - error: {1:e} - error/initial: {2:e} (|b-Ax|: {3:e})".format(i, error, error/1, self.get_grid_error())

            im.set_array(getarr())

            return im,
            #X = plt.imshow(getarr())
            #return X,

        ani = animation.FuncAnimation(fig, updatefig, interval=1, blit=False)
        plt.show()
        exit()

    def calc_grid_residual(self):
        self.grid.calc_residual(self.grid.w, self.grid.r, self.grid.u0, self.grid.u, matmul_func=self.grid.matmul)

    def get_grid_error(self):
        return self.grid.exactrro(self.grid.u0, self.grid.u, self.grid.matmul)

    def run_normal(self):
        cg_alphas = []
        cg_betas = []

        self.calc_grid_residual()
        initial_residual = np.sqrt(abs(self.get_grid_error()))

        # find E, initialise, etc
        self.grid.init_fadef(self.grid.u, self.grid.p, self.grid.r, self.grid.w, self.grid.z)

        #self.animate()

        rro = self.grid.dotvec(self.grid.r, self.grid.p)
        error = np.sqrt(abs(rro))

        i=0

        print "{0: 4d} - error: {1:e} - error/initial: {2:e} (|b-Ax|: {3:e})".format(i, error, error/initial_residual, self.get_grid_error())

        while error > self.params.tl_eps*initial_residual and i < self.params.tl_max_iters:
            if self.params.tl_use_cg:
                rro, alpha, beta = self.grid.cg(self.grid.u, self.grid.p, self.grid.r, self.grid.w, self.grid.z, rro, self.grid.matmul, M=self.grid.M)
            elif self.params.tl_use_dpcg:
                rro, alpha, beta = self.grid.fadef(self.grid.u, self.grid.p, self.grid.r, self.grid.w, self.grid.z)

            cg_alphas.append(alpha)
            cg_betas.append(beta)

            error = np.sqrt(abs(rro))

            i+=1
            if not i % 10:
                print "{0: 4d} - error: {1:e} - error/initial: {2:e} (|b-Ax|: {3:e})".format(i, error, error/initial_residual, self.get_grid_error())

        #print calc_eigs(cg_alphas, cg_betas)

        self.calc_grid_residual()
        error = np.sqrt(abs(self.get_grid_error()))

        print "EXACT error:", error/initial_residual, "at", i
        final_energy = self.grid.u[self.grid.inner]/self.grid.density[self.grid.inner]
        print "Energy:", np.linalg.norm(self.grid.energy)
        #self.plotu(self.grid.u)
        print

if __name__ == '__main__':
    solver = Solver()
    with Timer() as t:
        solver.run_normal()
    print t.interval


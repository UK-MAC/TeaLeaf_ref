import numpy as np
import matplotlib.pyplot as plt

def pltim(*arrs):
    #nplots = len(arrs)

    #plotscale=2
    #fig, axarr = plt.subplots(1, nplots)

    #if not isinstance(axarr, np.ndarray):
    #    axarr = [axarr]

    #for ii,(arr,ax) in enumerate(zip(arrs, axarr)):
    #    extent = [-1, arr.shape[0]+1, -1, arr.shape[1]+1]
    #    im = ax.imshow(arr, interpolation='nearest', cmap='spectral', extent=extent)
    #    plt.grid()
    #plt.colorbar(im)
    #plt.show()

    from mpl_toolkits.axes_grid1 import AxesGrid

    fig = plt.figure()

    grid = AxesGrid(fig, 111,
                    nrows_ncols = (1, len(arrs)),
                    axes_pad = 0.1,
                    label_mode = "1",
                    share_all = True,
                    cbar_location="bottom",
                    cbar_mode="each",
                    cbar_size="7%",
                    cbar_pad="2%",
                    )

    for i,arr in enumerate(arrs):
        extent = [0, arr.shape[0]+0, 0, arr.shape[1]+0]
        grid[i].title.set_text(i)
        im = grid[i].imshow(arr[1:-1, 1:-1], interpolation='nearest', cmap='seismic', extent=extent)
        grid.cbar_axes[i].colorbar(im)

    plt.show()

def plotu(*arrs):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter

    nplots = len(arrs)

    plotscale=2
    fig = plt.figure(figsize=(8*nplots*plotscale, 6*plotscale))

    for ii,arr in enumerate(arrs):
        uin = arr[1:-1, 1:-1]
        ax = fig.add_subplot(1, nplots, ii+1, projection='3d')
        X = np.arange(0, uin.shape[1])
        Y = np.arange(0, uin.shape[0])
        X, Y = np.meshgrid(X, Y)

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        cstride = uin.shape[0]/50
        rstride = uin.shape[1]/50

        shape = arr.shape
        ax.plot_surface(X, Y, uin, cmap=cm.coolwarm,
            cstride=cstride, rstride=rstride)

    plt.show()

class HasInner(object):
    # just to provide some convenience
    inner = np.s_[1:-1, 1:-1]

    xr_idx = np.s_[2:  , 1:-1]
    xl_idx = np.s_[ :-2, 1:-1]

    yr_idx = np.s_[1:-1, 2:  ]
    yl_idx = np.s_[1:-1,  :-2]

def calc_eigs(cg_alphas, cg_betas):
    cg_alphas = np.array(cg_alphas)
    cg_betas = np.array(cg_betas)

    offdiag = np.zeros_like(cg_alphas)
    diag = np.zeros_like(cg_alphas)

    diag = 1.0/cg_alphas

    diag[1:] += cg_betas[:-1]/cg_alphas[:-1]
    offdiag[1:] = np.sqrt(cg_betas[:-1])/cg_alphas[:-1]

    from scipy.linalg import eigvals_banded
    bands = np.vstack((offdiag, diag))
    eigs = eigvals_banded(bands)

    eigmin = eigs[0]
    eigmax = eigs[-1]

    #eigmin = eigmin*0.95
    #eigmax = eigmax*1.05

    return eigmin, eigmax

def calc_ch_coefs(eigmin, eigmax, bb, err, eps):
    cn = eigmax/eigmin
    print "Condition number:", cn

    #eps = 1e-10
    it_alpha = eps*bb/(4*err)
    gamm = (np.sqrt(cn) - 1)/(np.sqrt(cn + 1))
    est_itc = int(np.log(it_alpha)/(2*np.log(gamm)))
    print "Estimated itc:", est_itc

    theta = (eigmax + eigmin)/2
    delta = (eigmax - eigmin)/2

    ch_alphas = []
    ch_betas = []

    n_cheby_points = int(est_itc*10)

    sigma = theta/delta
    rho_old = 1.0/sigma
    for i in xrange(n_cheby_points):
        rho_new = 1.0/(2.0*sigma - rho_old)
        cur_alpha = rho_old*rho_new
        cur_beta = 2.0*rho_new/delta
        rho_old = rho_new

        ch_alphas.append(cur_alpha)
        ch_betas.append(cur_beta)

    ch_alphas = np.array(ch_alphas)
    ch_betas = np.array(ch_betas)

    return theta, ch_alphas, ch_betas


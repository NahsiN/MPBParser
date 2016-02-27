# plotting functions

from MPBParser import readfield, getscale
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.constants as spconsts
import numpy as np
import sys

def plotbands(mpb, bandlist=None, lw=1, xticks=None, xticklabels=None,
              figsize=None, ax_rect=None, has_light_line=False, ylims_offsets=[0, 0]):
    """
    Plots bands
    for light line assume ω = c|k|/n where n=1. In dimensionles coords ν = |k|
    """

    if figsize is not None:
        plt.figure(figsize=figsize)
    else:
        plt.figure()

    if ax_rect is not None:
        plt.axes(ax_rect)

    # Check if it makes sense to plot versus kmag
    if np.all(np.sort(mpb.kmag) == mpb.kmag):
        kindex_plot_flag = False
    else:
        kindex_plot_flag = True
        print('Nonsensical to use |k| for plotting.')

    # plot a specific number of bands
    if isinstance(bandlist, list):
        for band in bandlist:
            if kindex_plot_flag:
                plt.plot(mpb.freqs[:, band])
                if has_light_line:
                    plt.fill_between(range(len(mpb.kmag)), mpb.kmag, 1, alpha=0.5, facecolor='gray', edgecolor='black')
                    # plt.plot(range(len(mpb.kmag)), mpb.kmag)
            else:
                # plt.plot(mpb.kmag, mpb.freqs[:, band], '-b', linewidth=lw)
                plt.plot(mpb.kmag, mpb.freqs[:, band])
                if has_light_line:
                    plt.fill_between(mpb.kmag, mpb.kmag, 1, alpha=0.5, facecolor='gray', edgecolor='black')
                    # plt.plot(mpb.kmag, mpb.kmag)

        if kindex_plot_flag:
            plt.xlim(0, len(mpb.freqs[:, band])-1)
        else:
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])
            plt.xlabel(r'$|\mathbf{k}| \left[\frac{2\pi}{a}\right]$')

        # THIS 1e-3 OFFSET NEEDS TO BE TUNABLE FROM THE FUNCTION CALL
        plt.ylim(min(mpb.freqs[:, bandlist[-1]]) + ylims_offsets[0], max(mpb.freqs[:, bandlist[-1]]) + ylims_offsets[1])
        # plt.ylim(1e-3, )
        plt.ylabel(r'$\nu \left[\frac{c}{a}\right]$')
        # plt.tick_params(labelsize=ftsize)

    # plot all bands
    else:
        for band in range(mpb.numBands):
            if kindex_plot_flag:
                # plt.plot(mpb.freqs[:, band], color='b', linewidth=lw)
                plt.plot(mpb.freqs[:, band], color='b')
                ax = plt.gca()
                if xticklabels is not None:
                    ax.set_xticks(xticks)
                    ax.set_xticklabels(xticklabels)
                else:
                    print('You should specify x ticks manually.')
                    # ax.set_xticks((10, 20))
                    # ax.set_xticklabels(['Hi', 'Bye'])
                if has_light_line:
                    plt.fill_between(range(len(mpb.kmag)), mpb.kmag, 1, alpha=0.5, facecolor='gray', edgecolor='black')
                    # plt.plot(range(len(mpb.kmag)), mpb.kmag)

            else:
                # plt.plot(mpb.kmag, mpb.freqs[:, band], color='b', linewidth=lw)
                # plt.tick_params(labelsize=ftsize)
                plt.plot(mpb.kmag, mpb.freqs[:, band], color='b')
                if has_light_line:
                    plt.fill_between(mpb.kmag, mpb.kmag, 1, alpha=0.5, facecolor='gray', edgecolor='black')
                    # plt.plot(mpb.kmag, mpb.kmag)

        if kindex_plot_flag:
            plt.xlim(0, len(mpb.freqs[:, band])-1)
        else:
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])

        plt.ylim(1e-3, max(mpb.freqs[:, -1])+1e-2)
        plt.ylabel(r'$\nu \left[\frac{c}{a}\right]$')
        # plt.tick_params(labelsize=ftsize)


def plotfields(mpb, field_type, kindex=None, band=None, comp=None, mpbpostprocess=False,
    epsilon_contour_options={}, figsize=None):
    """
    Plots fields

    Inputs
    ------
    mpbpostprocess : False (default), assumes no post processing of the fields
                     has been performed using mpb-data
    """

    if figsize is not None:
        plt.figure(figsize=figsize)
    else:
        plt.figure()

    epsilon = readfield(mpb, field_type='epsilon_isotropic_trace', mpbpostprocess=mpbpostprocess)

    if field_type == 'e':
        E = readfield(mpb, kindex, band, field_type, mpbpostprocess=mpbpostprocess)
        E.create_complex()
        if comp is None:
            E2 = np.abs(E.x)**2 + np.abs(E.y)**2 + np.abs(E.z)**2
            if E2.ndim == 1:
                # (x) = getscale(mpb)
                # (xgrid, ygrid) = np.meshgrid(x, x)
                # creates a E2.ndim x E2.ndim square grid
                E2grid = np.tile(E2, (E2.shape[0], 1))
                epsilon_grid = np.tile(epsilon.dset, (epsilon.dset.shape[0], 1))
                # plt.contourf(xgrid, ygrid, E2grid)
                # plt.pcolormesh(E2grid)
                plt.imshow(E2grid)
                plt.colorbar()
                # plt.contour(xgrid, ygrid, epsilon_grid, colors='k', linewidths=lw)
                plt.contour(epsilon_grid, colors='w', **epsilon_contour_options)
                ax = plt.gca()
                ax.set_xticks(())
                ax.set_yticks(())
                plt.xlim(0, E2grid.shape[0]-1)
                plt.ylim(0, E2grid.shape[1]-1)

            elif E2.ndim == 2:
                # (x, y) = getscale(mpb)
                # (xgrid, ygrid) = np.meshgrid(x, y)
                # plt.contourf(xgrid, ygrid, E2)
                # plt.pcolormesh(E2)
                plt.imshow(E2)
                plt.colorbar()
                # plt.contour(xgrid, ygrid, epsilon.dset, colors='k', linewidths=lw)
                plt.contour(epsilon.dset, colors='w', **epsilon_contour_options)
                ax = plt.gca()
                ax.set_xticks(())
                ax.set_yticks(())
                # Determine which grid is on the axis. If Ny > Nx, the yindices
                # are placed on the x-axis
                if E2.shape[0] >= E2.shape[1]:
                    plt.xlim(0, E2.shape[0]-1)
                    plt.ylim(0, E2.shape[1]-1)
                else:
                    plt.xlim(0, E2.shape[1]-1)
                    plt.ylim(0, E2.shape[0]-1)

            elif E2.ndim == 3:
                # ASSUME SLAB GEOMETRY
                # xy cross section
                plt.imshow(E2[:, :, E2.shape[2]/2], aspect='equal')
                plt.colorbar()
                plt.contour(epsilon.dset[:, :, epsilon.dset.shape[2]/2], colors='w', **epsilon_contour_options)
                ax = plt.gca()
                ax.set_xticks(())
                ax.set_yticks(())

                # Determine which grid is on the axis. If Ny > Nx, the yindices
                # are placed on the x-axis
                if E2.shape[0] >= E2.shape[1]:
                    plt.xlim(0, E2.shape[0]-1)
                    plt.ylim(0, E2.shape[1]-1)
                else:
                    plt.xlim(0, E2.shape[1]-1)
                    plt.ylim(0, E2.shape[0]-1)

    elif field_type == 'epsilon':
        if len(epsilon.dset.shape) == 1:
            epsilon_grid = np.tile(epsilon.dset, (epsilon.dset.shape[0], 1))
            (x) = getscale(mpb)
            (xgrid, ygrid) = np.meshgrid(x, x)
            # plt.pcolormesh(xgrid, ygrid, epsilon_grid, cmap='Greys')
            plt.imshow(epsilon_grid, cmap='Greys')
            plt.colorbar()
            ax = plt.gca()
            ax.set_xticks(())
            ax.set_yticks(())
            plt.xlim(0, epsilon_grid.shape[0]-1)
            plt.ylim(0, epsilon_grid.shape[1]-1)

        elif len(epsilon.dset.shape) == 2:
            # Greys or gray
            # plt.pcolor(epsilon.dset, cmap='Greys')
            plt.imshow(epsilon.dset, cmap='Greys')
            plt.colorbar()
            plt.xlim(0, epsilon.dset.shape[0]-1)
            plt.ylim(0, epsilon.dset.shape[1]-1)
            ax = plt.gca()
            ax.set_xticks(())
            ax.set_yticks(())

        elif len(epsilon.dset.shape) == 3:
            # ASSUME PC SLAB GEOMETRY ONLY
            # plot in middle of slab
            # xy cross section
            plt.imshow(epsilon.dset[:, :, epsilon.dset.shape[2]/2], cmap='Greys', aspect='equal')
            plt.colorbar()
            # plt.xlim(0, epsilon.dset.shape[0]-1)
            # plt.ylim(0, epsilon.dset.shape[1]-1)
            ax = plt.gca()
            ax.set_xticks(())
            ax.set_yticks(())

    epsilon.close()


def plotvg(mpb, bandlist=None, lw=1, xticks=None, xticklabels=None,
              figsize=None, ax_rect=None, has_light_line=False):
    """
    Plots vg
    """

    if figsize is not None:
        plt.figure(figsize=figsize)
    else:
        plt.figure()

    if ax_rect is not None:
        plt.axes(ax_rect)

    # Check if it makes sense to plot versus kmag
    if np.all(np.sort(mpb.kmag) == mpb.kmag):
        kindex_plot_flag = False
    else:
        kindex_plot_flag = True
        print('Nonsensical to use |k| for plotting.')

    # plot a specific number of bands
    if isinstance(bandlist, list):
        for band in bandlist:
            # plt.plot(mpb.kmag, mpb.freqs[:, band], '-b', linewidth=lw)
            plt.plot(mpb.kmag, mpb.vgmag[:, band])

        if has_light_line:
            plt.fill_between(range(len(mpb.kmag)), mpb.kmag, 1, alpha=0.5, facecolor='gray', edgecolor='black')

        if kindex_plot_flag:
            plt.xlim(0, len(mpb.freqs[:, band])-1)
        else:
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])
            plt.xlabel(r'$|\mathbf{k}| \left[\frac{2\pi}{a}\right]$')
        plt.ylim(1e-3, np.max(mpb.vgmag) + 1e-2)
        # plt.ylim(1e-3, )
        plt.ylabel(r'$|v_g| [c]$')
        # plt.tick_params(labelsize=ftsize)

    # plot all bands
    else:
        for band in range(mpb.numBands):
            if kindex_plot_flag:
                # plt.plot(mpb.freqs[:, band], color='b', linewidth=lw)
                plt.plot(mpb.vgmag[:, band], color='b')
                ax = plt.gca()
                if xticklabels is not None:
                    ax.set_xticks(xticks)
                    ax.set_xticklabels(xticklabels)
                else:
                    print('You should specify x ticks manually.')
                    # ax.set_xticks((10, 20))
                    # ax.set_xticklabels(['Hi', 'Bye'])
            else:
                # plt.plot(mpb.kmag, mpb.freqs[:, band], color='b', linewidth=lw)
                # plt.tick_params(labelsize=ftsize)
                plt.plot(mpb.kmag, mpb.freqs[:, band], color='b')
        if has_light_line:
            plt.fill_between(range(len(mpb.kmag)), mpb.kmag, 1, alpha=0.5, facecolor='gray', edgecolor='black')

        if kindex_plot_flag:
            plt.xlim(0, len(mpb.vgmag[:, band])-1)
        else:
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])
            plt.xlabel(r'$|\mathbf{k}| \left[\frac{2\pi}{a}\right]$')

        plt.ylim(1e-3, np.max(mpb.vgmag)+1e-2)
        plt.ylabel(r'$|v_g| \left[c\right]$')
        # plt.tick_params(labelsize=ftsize)

# -*- coding: utf-8 -*-
# plotting functions

from MPBParser import readfield, getscale
from scipy.fftpack import fftshift, fft2
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.constants as spconsts
import numpy as np
import sys
import ipdb

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
                ax = plt.gca()
                if xticklabels is not None:
                    ax.set_xticks(xticks)
                    ax.set_xticklabels(xticklabels)
                else:
                    print('You should specify x ticks manually.')
                if has_light_line:
                    plt.fill_between(range(len(mpb.kmag)), mpb.kmag, 1, alpha=1.0, facecolor='gray', edgecolor='black', zorder=3)
                    # plt.plot(range(len(mpb.kmag)), mpb.kmag)
            else:
                # plt.plot(mpb.kmag, mpb.freqs[:, band], '-b', linewidth=lw)
                plt.plot(mpb.kmag, mpb.freqs[:, band])
                if has_light_line:
                    plt.fill_between(mpb.kmag, mpb.kmag, 1, alpha=1.0, facecolor='gray', edgecolor='black', zorder=3)
                    # plt.plot(mpb.kmag, mpb.kmag)

        if kindex_plot_flag:
            plt.xlim(0, len(mpb.freqs[:, band])-1)
        else:
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])
            plt.xlabel(r'$|\mathbf{k}| \left[\frac{2\pi}{a}\right]$')

        # THIS 1e-3 OFFSET NEEDS TO BE TUNABLE FROM THE FUNCTION CALL
        plt.ylim(np.min(mpb.freqs[:, bandlist]) + ylims_offsets[0], np.max(mpb.freqs[:, bandlist]) + ylims_offsets[1])
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
                    plt.fill_between(range(len(mpb.kmag)), mpb.kmag, 1, alpha=1.0, facecolor='gray', edgecolor='black', zorder=3)
                    # plt.plot(range(len(mpb.kmag)), mpb.kmag)

            else:
                # plt.plot(mpb.kmag, mpb.freqs[:, band], color='b', linewidth=lw)
                # plt.tick_params(labelsize=ftsize)
                plt.plot(mpb.kmag, mpb.freqs[:, band], color='b')
                if has_light_line:
                    plt.fill_between(mpb.kmag, mpb.kmag, 1, alpha=1.0, facecolor='gray', edgecolor='black', zorder=3)
                    # plt.plot(mpb.kmag, mpb.kmag)

        if kindex_plot_flag:
            plt.xlim(0, len(mpb.freqs[:, band])-1)
        else:
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])
            plt.xlabel(r'$|\mathbf{k}| \left[\frac{2\pi}{a}\right]$')

        plt.ylim(ylims_offsets[0], max(mpb.freqs[:, -1]) + ylims_offsets[1])
        plt.ylabel(r'$\nu \left[\frac{c}{a}\right]$')
        # plt.tick_params(labelsize=ftsize)


def plotfields(mpb, field_type, kindex=None, band=None, comp=None, mpbpostprocess=False,
    epsilon_contour_options={}, figsize=None, field_file=None, epsilon_file=None, plot_options={}):
    """
    Plots fields

    Inputs
    ------
    mpbpostprocess : False (default), assumes no post processing of the fields
                     has been performed using mpb-data
    field_type : 'e', 'epsilon', 'epsilon-ft'
    plot_options : specifies extra plot options such as axes limits
    comp: 'x', 'y', 'z'. Only rectangular coordinates for now
    """

    if figsize is not None:
        plt.figure(figsize=figsize)
    else:
        plt.figure()

    if epsilon_file is None:
        epsilon = readfield(mpb, field_type='epsilon_isotropic_trace', mpbpostprocess=mpbpostprocess)
    elif isinstance(epsilon_file, str):
        epsilon = readfield(mpb, field_type='epsilon_isotropic_trace', mpbpostprocess=mpbpostprocess, field_file=epsilon_file)

    if field_type == 'e':
        if field_file is None:
            E = readfield(mpb, kindex, band, field_type, mpbpostprocess=mpbpostprocess)
        elif isinstance(field_file, str):
            E = readfield(mpb, kindex, band, field_type, mpbpostprocess=mpbpostprocess, field_file=field_file)

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
                plt.colorbar(label=r'$|\mathbf{E}|^2$')
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
                plt.colorbar(label=r'$|\mathbf{E}|^2$')
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
                plt.colorbar(label=r'$|\mathbf{E}|^2$')
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
        else:
            if comp == 'x':
                E_comp2 = np.abs(E.x)**2
            elif comp == 'y':
                E_comp2 = np.abs(E.y)**2
            elif comp == 'z':
                E_comp2 = np.abs(E.z)**2
            else:
                raise('Invalid component entered. x,y,z are valid options')

            if E_comp2.ndim == 1:
                pass
            elif E_comp2.ndim == 2:
                pass
            elif E_comp2.ndim == 3:
                # ASSUME SLAB GEOMETRY
                # xy cross section
                plt.imshow(E_comp2[:, :, E_comp2.shape[2]/2], aspect='equal')
                plt.colorbar(label=r'$|\mathbf{E}|^2$')
                plt.contour(epsilon.dset[:, :, epsilon.dset.shape[2]/2], colors='w', **epsilon_contour_options)
                ax = plt.gca()
                ax.set_xticks(())
                ax.set_yticks(())

                # Determine which grid is on the axis. If Ny > Nx, the yindices
                # are placed on the x-axis
                if E_comp2.shape[0] >= E_comp2.shape[1]:
                    plt.xlim(0, E_comp2.shape[0]-1)
                    plt.ylim(0, E_comp2.shape[1]-1)
                else:
                    plt.xlim(0, E_comp2.shape[1]-1)
                    plt.ylim(0, E_comp2.shape[0]-1)



    elif field_type == 'epsilon':
        if len(epsilon.dset.shape) == 1:
            epsilon_grid = np.tile(epsilon.dset, (epsilon.dset.shape[0], 1))
            # epsilon_grid_ft = fftshift(fft2(epsilon_grid))
            (x) = getscale(mpb)
            (xgrid, ygrid) = np.meshgrid(x, x)
            # plt.pcolormesh(xgrid, ygrid, epsilon_grid, cmap='Greys')
            plt.imshow(epsilon_grid, cmap='Greys')
            plt.colorbar(label=r'$\varepsilon$')
            ax = plt.gca()
            ax.set_xticks(())
            ax.set_yticks(())
            plt.xlim(0, epsilon_grid.shape[0]-1)
            plt.ylim(0, epsilon_grid.shape[1]-1)
            # plt.figure()
            # plt.imshow(np.abs(epsilon_grid_ft)**2)
            # plt.colorbar()

        elif len(epsilon.dset.shape) == 2:
            # Greys or gray
            # plt.pcolor(epsilon.dset, cmap='Greys')
            # filter out DC_component
            # epsilon_ft[np.unravel_index(np.argmax(epsilon_ft), np.shape(epsilon_ft))] = 0

            plt.imshow(epsilon.dset, cmap='Greys')
            plt.colorbar(label=r'$\varepsilon$')
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
            # yz cross section
            # plt.imshow(epsilon.dset[epsilon.dset.shape[0]/2, :, :], cmap='Greys', aspect='equal')
            # plt.imshow(epsilon.dset[0, :, :], cmap='Greys', aspect='equal')
            plt.colorbar(label=r'$\varepsilon$')
            # plt.xlim(0, epsilon.dset.shape[0]-1)
            # plt.ylim(0, epsilon.dset.shape[1]-1)
            ax = plt.gca()
            ax.set_xticks(())
            ax.set_yticks(())

    elif field_type == 'epsilon-ft':
        if len(epsilon.dset.shape) == 1:
            pass
        elif len(epsilon.dset.shape) == 2:
            # zero pad
            epsilon_ft = fftshift(fft2(np.pad(epsilon.dset[:,:], (512, 512), 'constant')))
            plt.imshow((np.abs(epsilon_ft)/np.max(np.abs(epsilon_ft)))**2, cmap='Greys', vmin=0, vmax=0.1)
            ax = plt.gca()
            ax.set_xticks(())
            ax.set_yticks(())

        elif len(epsilon.dset.shape) == 3:
            pass

        ax = plt.gca()
        ax.set_xticks(())
        ax.set_yticks(())

        if 'xlims' in plot_options:
            xa = plot_options['xlims'][0]
            xb = plot_options['xlims'][1]
            plt.xlim(xa, xb)

        if 'ylims' in plot_options:
            ya = plot_options['ylims'][0]
            yb = plot_options['ylims'][1]
            plt.ylim(ya, yb)

    epsilon.close()


def plotvg(mpb, bandlist=None, lw=1, xticks=None, xticklabels=None,
              figsize=None, ax_rect=None, has_light_line=False, ylims_offsets=[0, 0]):
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
            plt.fill_between(range(len(mpb.kmag)), mpb.kmag, 1, alpha=1.0, facecolor='gray', edgecolor='black')

        if kindex_plot_flag:
            plt.xlim(0, len(mpb.freqs[:, band])-1)
        else:
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])
            plt.xlabel(r'$|\mathbf{k}| \left[\frac{2\pi}{a}\right]$')

        # plt.ylim(1e-3, np.max(mpb.vgmag) + 1e-2)
        plt.ylim(np.min(mpb.vgmag[:, bandlist]) + ylims_offsets[0], np.max(mpb.vgmag[:, bandlist]) + ylims_offsets[1])
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
            plt.fill_between(range(len(mpb.kmag)), mpb.kmag, 1, alpha=1.0, facecolor='gray', edgecolor='black')

        if kindex_plot_flag:
            plt.xlim(0, len(mpb.vgmag[:, band])-1)
        else:
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])
            plt.xlabel(r'$|\mathbf{k}| \left[\frac{2\pi}{a}\right]$')

        plt.ylim(ylims_offsets[0], np.max(mpb.vgmag) + ylims_offsets[1])
        plt.ylabel(r'$|v_g| \left[c\right]$')
        # plt.tick_params(labelsize=ftsize)

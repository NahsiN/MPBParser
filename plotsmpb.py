# plotting functions

import matplotlib.pyplot as plt
import scipy.constants as spconsts

plt.rc('text', usetex=True)


def plotbands(mpb, bandnum=None, ftsize=10, lw=1):
    """
    Plots bands
    """

    plt.figure()
    if isinstance(bandnum, int):
        plt.plot(mpb.kmag, mpb.freqs[:, bandnum], linewidth=lw)
        plt.xlim(mpb.kmag[0], mpb.kmag[-1])
        # plt.ylim(1e-3, )
        plt.ylabel(r'$\nu [\frac{c}{a}]$', fontsize=ftsize)
        plt.tick_params(labelsize=ftsize)
    else:
        for band in range(mpb.numBands):
            plt.plot(mpb.kmag, mpb.freqs[:, band], color='b', linewidth=lw)
            plt.xlim(mpb.kmag[0], mpb.kmag[-1])
            plt.ylabel(r'$\nu [\frac{c}{a}]$', fontsize=ftsize)
            plt.tick_params(labelsize=ftsize)
            plt.ylim(1e-3, max(mpb.freqs[:, -1])+1e-3)


def plotfields():
    """
    Plots fields
    """
    
    pass

# -*- coding: utf-8 -*-
import numpy as np
import math
import scipy.constants as spconsts
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq, fftshift
from scipy import signal
import sys
import ipdb
import h5py

I = complex(0, 1)

def toSI(x, quantity_type, SI_units, norm_length=1, norm_time=1, norm_speed=1):
    """
    Convert normalized units to SI units.

    Input
    -----
    x : a normalized quantity, a number, a numpy array etc. list or tuple
    quantity_type : 'freq', 'length', 'time'
    SI_units : units of desired quantity.
               Frequency : ('Hz', 'GHz', 'THz')
               Length : ('m', 'mm', 'microm', 'nm')
               Time : ('fs')
    norm_length : normalization units for length
    norm_time : normalization units for time
    norm_speed : normalization units for speed

    Returns
    -------
    (x_SI, 'units') : tuple containing the quantity/quantities in SI and it's corresponding unit

    Examples
    --------
    (xSI, units) = toSI(x, 'freq', 'thz', norm_length=L, norm_time=T, norm_speed=L/T)
    x_SI = toSI(x, 'freq', 'thz', norm_length=L, norm_time=T, norm_speed=L/T) where
    x_SI[0] is the number of x_SI[1] is the unit
    """
    norm_units = set([norm_speed, norm_length, norm_time])  # set contains no duplicate elements

    c = norm_speed
    L = norm_length
    T = norm_time

    x = np.array(x)  # cast x as a numpy array to allow for element-wise operations
    if len(norm_units) == 1:
        print('No conversion of units required. All fundamental units are same')
        return None
    elif c != L/T:
        print('c≠L/T. Look at your normalized units')
        return None

    if quantity_type == 'freq':
        #x_SI = [x*c/L, 'Hz']
        if SI_units == 'THz' or SI_units == 'thz':
            x_SI = (x*c/L/spconsts.tera, 'THz')
        elif SI_units == 'GHz' or SI_units == 'ghz':
            x_SI = (x*c/L/spconsts.giga, 'GHz')
        else:
            x_SI = (x*c/L, 'Hz')

    elif quantity_type == 'length':
        if SI_units == 'nm':
            x_SI = (x*L/spconsts.nano, 'nm')
        elif SI_units == 'microm':
            x_SI = (x*L/spconsts.micro, 'μm')
        elif SI_units == 'mm':
            x_SI = (x*L/spconsts.milli, 'mm')
        else:
            x_SI = (x*L, 'm')

    elif quantity_type == 'time':
        if SI_units == 'fs':
            x_SI = (x*T/spconsts.femto, 'fs')
        elif SI_units == 'ps':
            x_SI = (x*T/spconsts.pico, 'ps')
        else:
            x_SI = (x*T, 's')

    else:
        print('The quantity type entered is not valid, valid options are freq, length, time')
        return None

    return x_SI

def nextpow2(N, npow=0):
    """
    Returns integer N' s.t N'=2^x where x is an integer
    """
    N_prime = np.ceil(np.log2(N))  # N_prime is a power of two
    return(int(2**(N_prime + npow)))


def fourier_transform(funcs, dts, inSI=False, SI_units='thz',
                    convention='math', zero_padding=False, npow=1,
                     **kwargs_norm_units):
    """
    The F.T. of an array with various sampling times for each of the funcs

    Input
    -----
    funcs : a list of np.array types representing the discrete time signal.
    I assume you have correctly sampled the funcs by not including the endpoint
    of the time signal
    Must be casted as list so enclose with []. Or as a tuple like (f1,) with the comma being important. CAST AS A TUPLE
    dts : the sampling time step. For now, every function has the same time step.
    inSI (optional) : wheteher to use SI units or normalized units for frequency (Default: False)
    SI_units (optional) : specify the SI unit for frequency. (Default: THz). Valid options are 'thz',
    **kwargs_norm_units (optional) : If inSI==True, then specify the appropriate normalized units for
    the toSI function. Please see docstring for toSI for further details.
    convention: 'physics' convention or 'math' convention. Physics convention uses
    ifft as Fourer transform and math convention uses fft
    Returns
    -------
    (frequencies, fft of signal) :  list containing the frequencies and
    spectrum with the zero frequency centered
    else (Default)
    Return nothing
    """

    if isinstance(dts, float):
        dts = np.tile(dts, len(funcs))  # cast as a list
    elif isinstance(dts, list):
        if len(dts) != len(funcs):
            print('The funcs and dts array are not the same size.')
            return None

    funcs_freq = []
    freqs_list = []
    for (f, dt) in zip(funcs, dts):
        Ns = np.size(f)
        # windowing
        # window_real = np.max(np.real(f))*signal.hann(Ns)
        # window_imag = np.max(np.imag(f))*signal.hann(Ns)
        # Testing window
        # plt.figure()
        # plt.plot(window)
        # plt.plot(np.abs(f))
        # plt.show()
        # plt.close()
        # f = window*f
        # f = window_real*np.real(f) + I*window_imag*np.imag(f)
        if zero_padding:
            NFFT = nextpow2(Ns, npow=npow)
            # pad asymmetrical on right side, assumes sampling in [0,T)
            # f = np.append(f, np.zeros(NFFT-Ns))
            # f = np.pad(f, (0, NFFT-Ns), 'constant', constant_values=(0, 0))

            # pad symmetrically on both sides since it's a gaussian. sampling [-0.5*T, 0.5*T)
            Nleft = int(np.ceil(0.5*(NFFT - Ns)))
            Nright = int((NFFT - Ns) - Nleft)
            f = np.pad(f, (Nleft, Nright), 'constant', constant_values=(0, 0))
            Ns = np.size(f)
        # ipdb.set_trace()
        if convention == 'math':
            f_freq = fftshift((fft(f)))
        elif convention == 'physics':
            # f_freq = fftshift((ifft(f, n=2**15)))  # another exaple of asymmetric padding on the right side
            # Ns = np.size(f_freq)
            f_freq = fftshift((ifft(f)))
        freqs = fftshift(fftfreq(Ns, dt))

        funcs_freq.append(f_freq)
        freqs_list.append(freqs)
        if inSI == True:
            freqs = toSI(freqs, 'freq', SI_units, **kwargs_norm_units)
            if freqs == None:
               print('toSI function call failed. Aborting')
               return None
    return (freqs_list, funcs_freq)


def find_indices(x, values, tol=1e-3):
    """
    Takes in an array and a list of values and tol. Returns the indices
    Inputs
    ------
    x : array of values
    values : values we are looking for in x
    tol : accepted tolerance

    Returns
    -------
    indices : a list of indices containing the values
    """

    indices = []
    for xi in values:
        # Find indices which closely match the x points you need
        # check if x is hdf dataset
        if isinstance(x, h5py._hl.dataset.Dataset):
            xinds = np.argwhere(abs(x - xi) <= tol)  # find all the
        else:
            xinds = np.argwhere(abs(x - xi) <= tol)  # find all the non-zeros (True) xinds

        if not len(xinds) > 0:
            print('Please decrease your tol={0:.3e} or you have entered a x point which is not in the x array'.format(tol))
            print('No matches to x_ideal={0:.3f} were found'.format(xi))
            sys.exit()
        xind = xinds[0, 0]  # just select the first one. if only [0] is used one gets an array. np.argwhere outputs array of arrays
        indices.append(xind)
    return indices


def latexfigsize(width, columnwidth_latex=433.62):
    """
    Reference: http://scipy-cookbook.readthedocs.io/items/Matplotlib_LaTeX_Examples.html#producing-graphs-for-publication-using-latex
    """
    # Get this from LaTeX using \showthe\columnwidth
    fig_width_pt = width*columnwidth_latex
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    golden_mean = (math.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    # fig_height = fig_width*golden_mean       # height in inches
    fig_height = 0.5*golden_mean*fig_width     # for 1x2 subplots
    fig_size = (fig_width, fig_height)
    return fig_size

def freq_filter(f, dt, cutoff_freq, convention='math'):
    """
    A digital filter that filters frequency below a cutoff frequency

    Parameters
    ----------
    f : time signal
    dt : sampling period
    nu_cutoff : cutoff frequency

    Returns
    -------
    The filtered time signal
    """

    if convention == 'math':
        f_freq = fft(f)
    elif convention == 'physics':
        f_freq = ifft(f)
    Ns = np.size(f)
    freqs = fftfreq(Ns, dt)

    # filtering operation
    f_freq[np.where(np.abs(freqs) > cutoff_freq)] = 0
    # go back to time domain
    if convention == 'math':
        f_filter_time = ifft(f_freq)
    elif convention == 'physics':
        f_filter_time = fft(f_freq)

    return f_filter_time

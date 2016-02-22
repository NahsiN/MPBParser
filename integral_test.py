from MPBParser import MPBBandStructure, readfield, getscale
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spintegrate
from indicator import indicator_func

mpb = MPBBandStructure('/home/nishan/Code/thales/MPB/w14/nonbloch/w14.out', symm='zeven')
# mpb.csvparser()
mpb.readbanddata()
#E_nonbloch = readfield(mpbobj=mpb, field_type='e.nonbloch.v', kindex=3, band=16)
E_nonbloch = readfield(field_file='/home/nishan/Code/thales/MPB/w14/nonbloch/e.nonbloch.v.k10.b16.zeven.h5', field_type='e.nonbloch.v')
E_nonbloch.create_complex()

# (x, dx, y, dy, z, dz) = getscale(mpb, retstep=True)

ind = indicator_func(mpb, 10.0489, 1.0, 'slab_only')
# for DEBUGGING
# plt.figure()
# plt.imshow(ind[:,:, 50])
# plt.colorbar()
# plt.show()


def nlo_coeffs(efield, indicator, ccnl):
    """
    Returns the nonlinear coupling coeffs after integration
    efield is assumed to be nonbloch
    """

    (x, y, z) = getscale(mpb)
    # (x, dx, y, dy, z, dz) = getscale(mpb, retstep=True)
    Nx = np.size(x)
    # Ny = np.size(y)
    # Nz = np.size(z)
    # Symbols : ((+,+), |+|), ((+,+), |-|)
    # compute the field quantities
    e2 = np.abs(efield.x)**2 + np.abs(efield.y)**2 + np.abs(efield.z)**2
    e4 = e2**2
    ex4 = np.abs(efield.x)**4
    ey4 = np.abs(efield.y)**4
    ez4 = np.abs(efield.z)**4
    e2complex = efield.x**2 + efield.y**2 + efield.z**2
    ex2 = np.abs(efield.x)**2
    ey2 = np.abs(efield.y)**2
    ez2 = np.abs(efield.z)**2
    ex2complex = efield.x**2
    ey2complex = efield.y**2
    ez2complex = efield.z**2

    for cc_type in ccnl.keys():
        cc = np.zeros(Nx, dtype=complex)
        print('Solving for NL coefficient {0}'.format(cc_type))
        if cc_type == (('+', '+'), '|+|'):
            integrand = indicator*(e4 + 2*(ex4 + ey4 + ez4))
        elif cc_type == (('+', '+'), '|-|'):
            integrand = indicator*e4
        elif cc_type == (('+', '-'), ('+', '-')):
            integrand = indicator*(ex4 + ey4 + ez4)
        elif cc_type == (('+', '+'), ('+', '-')):
            integrand = indicator*((e2complex + 2*ex2complex)*ex2 + (e2complex + 2*ey2complex)*ey2 + (e2complex + 2*ez2complex)*ez2)
        elif cc_type == (('+', '-'), '|-|'):
            integrand = indicator*((e2 + 2*ex2)*np.conjugate(ex2complex) + (e2 + 2*ey2)*np.conjugate(ey2complex) + (e2 + 2*ez2)*np.conjugate(ez2complex))
        elif cc_type == (('+', '-'), '|+|'):
            integrand = indicator*(e2*np.conjugate(ex2complex + ey2complex + ez2complex))
        elif cc_type == (('+', '+'), ('-', '+')):
            integrand = indicator*(np.conjugate(ex2complex)*ex2 + np.conjugate(ey2complex)*ey2 + np.conjugate(ez2complex)*ez2)
        elif cc_type == (('+', '-'), ('-', '+')):
            integrand = indicator*(np.conjugate((e2complex + 2*ex2complex)*ex2complex + (e2complex + 2*ey2complex)*ey2complex + (e2complex + 2*ez2complex)*ez2complex))

        for i in range(Nx):
            cc[i] = spintegrate.simps(spintegrate.simps(integrand[i,:,:], z, axis=1), y, axis=0)
            # cc[i] = spintegrate.simps(spintegrate.simps(integrand[i,:,:], dx=dz, axis=1), dx=dy, axis=0)

        ccnl[cc_type] = cc  # modified in place. [:] is NECESSARY

    print('Done!')

# cc1 = nlo_coeffs(E_nonbloch, indicator=ind, symbol='((+,+),|+|)')
ccnl = {(('+', '+'), '|+|'): None, (('+', '+'), '|-|'): None, (('+', '-'), ('+', '-')): None,
        (('+', '+'), ('+', '-')): None, (('+', '-'), '|-|'): None, (('+', '-'), '|+|'): None,
        (('+', '+'), ('-', '+')): None, (('+', '-'), ('-', '+')): None}
nlo_coeffs(E_nonbloch, ind, ccnl)
# cc2 = nlo_coeffs(E_nonbloch, indicator=ind, symbol='((+,+),|-|)')

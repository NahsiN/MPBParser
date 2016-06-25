# test hole integral step function
from MPBParser import MPBBandStructure, readfield
import numpy as np
import matplotlib.pyplot as plt

#mpb = MPBBandStructure('/home/nishan/Code/thales/MPB/w14/w14.out', 'zeven')
#mpb.csvparser()
#mpb.readbanddata()

#eps_slab = 10.0489
#eps_air = 1.0
#TOL = 1e-3
# no need since I start with an array of zeros
# mask_air_holes_supercell = abs(epsilon.dset[:]- eps_air) <= TOL
# step_func[mask_air_holes_supercell] = 0

def indicator_func(mpb, eps_slab, eps_air, type, TOL=1e-3):
    # type air_slab_bdry, slab_only, air_hole_slab_bdry
    epsilon = readfield(mpb, field_type='epsilon_isotropic_trace')

    if type == 'slab_only':
        indicator_slab = np.zeros(epsilon.dset.shape)
        mask_slab = abs(epsilon.dset[:] - eps_slab) <= TOL
        indicator_slab[mask_slab] = 1
        indicator = indicator_slab


    elif type == 'air_slab_bdry':
        indicator_air_slab_bdry = np.zeros(epsilon.dset.shape)
        # still keeps the slab-hole boundary
        mask_air_slab_bdry = np.logical_not(np.logical_or(abs(epsilon.dset[:]- eps_air) <= TOL, \
            abs(epsilon.dset[:]- eps_slab) <= TOL))
        indicator_air_slab_bdry[mask_slab_hole_bdry] = 1
        indicator = indicator_air_slab_bdry

    elif type == 'air_hole_slab_bdry':
        indicator_air_hole_slab_bdry = np.zeros(epsilon.dset.shape)
        # still keeps the slab-hole boundary
        mask_air_hole_slab_bdry = np.logical_not(np.logical_or(abs(epsilon.dset[:]- eps_air) <= TOL, \
            abs(epsilon.dset[:]- eps_slab) <= TOL))
        # loop over the z dimension and weed out slab-supercell boundary by looking
        # at mean epsilon in the xy cross section
        for k in range(epsilon.dset.shape[2]):
            # 0.5*(eps_air + eps_slab) is a reasonable estimate.
            if np.mean(epsilon.dset[:,:,k]) < 0.5*(eps_air + eps_slab):
                mask_air_hole_slab_bdry[:,:,k] = 0
        indicator_air_hole_slab_bdry[mask_slab_hole_bdry] = 1
        indicator = indicator_air_hole_slab_bdry
        # FOR DEBUGGING
        # plt.figure()
        # plt.imshow(indicator_slab_hole_bdry[:,:, int(epsilon.dset.shape[2]/2)])
        # plt.colorbar()
        # plt.show()
    else:
        print('Invalid type entereted. Valid options are')
        epsilon.close()
        return None
    epsilon.close()
    return indicator

# reads an .h5 file from disk

import h5py
import numpy as np

fname = '/home/nishan/Code/thales/MPB/w14/e.k03.b16.zeven.h5'
f = h5py.File(fname, 'r')

# ex_real = f['x.r']

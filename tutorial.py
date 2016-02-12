from MPBParser import MPBBandStructure, readfield

mpb = MPBBandStructure('/home/nishan/Code/thales/MPB/w14/w14.out', 'zeven')
mpb.csvparser()
mpb.readbanddata()
E = readfield(mpb,field_type='e',kindex=3,band=16)
H = readfield(mpb,field_type='h',kindex=3,band=16)
epsilon_tensor = readfield(mpb, field_type = 'epsilon')
epsilon = readfield(mpb, field_type = 'epsilon_isotropic_trace')

# CLOSE THE H5 FILE OBJECTS WHEN DONE!
# E_h5_fobj.close()

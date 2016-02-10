from MPBParser import MPBBandStructure, readfield

mpb = MPBBandStructure('/home/nishan/Code/thales/MPB/w14/w14.out', 'zeven')
mpb.csvparser()
mpb.readbanddata()
E = readfield(3,16,mpb)
H = readfield(3,16,mpb,field_type='h')

# CLOSE THE H5 FILE OBJECTS WHEN DONE!
# E_h5_fobj.close()

## Numerical things
- [ ] **Med** stop the overeliance on mpb object
- [ ] **Low** tree structure for mpb geomteric objects
- [ ] **Low** get ε for slab and air using tree structure

## Error Checks
- [ ] **Low** in readfield error for invalid field_type
- [ ] **Medium** Check for invalid filename in readfield
- [ ] **Low** getscale rectangular grid check
- [ ] **Low** MPBParser won't work with -k

## Completed
- [x] lightline included if needed
- [x] .csv file not created if present
- [x] **Low** plot the bands, fields, epsilon
- [x] mpbpostprocess flag allows plotting after running of mpb-data
- [x] 2D plot of field and bandstructure
- [x] checks to see if it is sensible to plot against kmag
- [x] .csv files are created in mpb.path
- [x] **High** If MPB out doesn't have velocities, don't bother creating data attrs related to velocity
- [x] **Low** Plot bands basic
- [x] **Med** strip phase from EM fields if needed. MPB does it!
- [x] **High** create grid. Equivalent to getScale
- [x] **High** do an integral
- [x] **Low** readbandata can determine numK, numBands from csv file alone
- [x] MPBBandstructure class to load all the band data from mpb output
- [x] readfield function to read fields from mpb output
- [x] **High** indicator function for slab only, slab-hole boundary only
- [x] **Med** base class for tensor types having the close() function in common

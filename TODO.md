## Numerical things
- [ ] **Med** stop the overeliance on mpb object
- [ ] **Low** plot the bands, fields, epsilon
- [ ] **Low** tree structure for mpb geomteric objects
- [ ] **Low** get ε for slab and air using tree structure
- [ ] **High** If MPB out doesn't have velocities, don't bother

## Error Checks
- [ ] **Low** in readfield error for invalid field_type
- [ ] **Medium** Check for invalid filename in readfield
- [ ] **Low** getscale rectangular grid check


## Completed
- [x] **Low** Plot bands basic
- [x] **Med** strip phase from EM fields if needed. MPB does it!
- [x] **High** create grid. Equivalent to getScale
- [x] **High** do an integral
- [x] **Low** readbandata can determine numK, numBands from csv file alone
- [x] MPBBandstructure class to load all the band data from mpb output
- [x] readfield function to read fields from mpb output
- [x] **High** indicator function for slab only, slab-hole boundary only
- [x] **Med** base class for tensor types having the close() function in common

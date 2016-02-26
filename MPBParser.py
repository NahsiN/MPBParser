# Parse the MPB output file for relevant quantities
# Inspired from Patterson's Matlab code

# parse the text file
import re
import csv
import numpy as np
import h5py
import sys
import os

def strip_spaces(fields):
    ''' Strip spaces and newline characters from a list of strings '''
    # id(fields[:]) ≠ id(fields). fields[:] creates a new copy to allow
    # me to edit fields inside the loop (Remember Python passes by reference)
    # strip out spaces
    for element in fields[:]:
        fields.pop(0)
        fields.append(element.strip(' \n'))  # \n only needed for last element
    return fields


class MPBBandStructure:
    """
    Class for parsing MPB output
    """

    def __init__(self, out_file, symm):
        self.fname_mpb = out_file
        self.symm = symm
        # self.path = out_file.rstrip(re.search('[\w.]+.out', out_file).group())  # relies on file being named foo.out
        self.path = '/'.join(out_file.split('/')[0:-1])  # [0:-1] ensures not to include filename. Spaces allowed?

    def csvparser(self):
        """
        Parses MPB output file to two csv files while defining class data attributes
        """

        symm = self.symm
        fname_mpb = self.fname_mpb

        f_mpb = open(fname_mpb, 'r')
        f_freqs_path = '/'.join([self.path, symm+'freqs.csv'])
        f_velocities_path = '/'.join([self.path, self.symm+'velocities.csv'])

        # Check to see if a .csv file already exists
        if os.path.isfile(f_freqs_path) and os.path.isfile(f_velocities_path):
            print('{0} and {1} already exist. SKIP CREATION.'.format(f_freqs_path, f_velocities_path))
            return None

        print('Creating {0}'.format(f_freqs_path))
        print('Creating {0}'.format(f_velocities_path))
        f_freqs = open(f_freqs_path, 'w', newline='')
        f_velocities = open(f_velocities_path, 'w', newline='')

        freqs_writer = csv.writer(f_freqs)  # initialize a writer object
        velocities_writer = csv.writer(f_velocities)

        for line in f_mpb:
            # print(line)  # DEBUG
            # determine number of k points and bands
            if re.match(r'Computing [\d]* bands', line):
                self.numBands = int(re.findall(r'[\d]+', line)[0])
            if re.match(r'[\d]+ k-points', line):
                self.numK = int(re.findall(r'[\d]+', line)[0])

            # determine a header to write to the velocites .csv file
            # search for the line 'zevenfreqs:, k index, k1, .....'
            if re.search(symm + r'freqs:, k index', line):
                # write velocites csv file header
                header = []
                header.append('k index')
                header_rest = (re.search(symm+r' band [\s\S]*', line).group()).split(',')
                for element in header_rest:
                    header.append(element)
                header = strip_spaces(header)
                velocities_writer.writerow(header)  # write out list as row in .csv file

            # look for lines of the form 'zevenfreqs:,'
            # write out to frequencies .csv file
            if re.match(symm + r'freqs:,[\s\S]*', line):
                # create a list by splitting line
                fields = (line.lstrip(symm+'freqs:,')).split(',')
                fields = strip_spaces(fields)
                freqs_writer.writerow(fields)  # write out list as row in .csv file

            # write out velocities to .csv file
            if re.match(symm + r'velocity:,[\s\S]*', line):
                fields = (line.lstrip(symm+'velocity:,')).split(',')
                for element in fields[:]:
                    fields.pop(0)
                    fields.append(element.strip(' #()\n'))
                velocities_writer.writerow(fields)  # write out list as row in .csv file

        f_mpb.close()
        f_freqs.close()
        f_velocities.close()

    def readbanddata(self):
        """
        Loads csv files to parse their textual data into numpy arrays. Creates
        new data attributes for object.

        self.k: [self.numK × 3]
        self.kmag: [self.numK × 1]
        self.freqs: [self.numK × self.numBands]
        self.vg: [self.numK × 3 × self.numBands]
        self.vgmag: [self.numK × self.numBands]
        """

        symm = self.symm

        f_freqs_path = '/'.join([self.path, symm+'freqs.csv'])
        f_velocities_path = '/'.join([self.path, self.symm+'velocities.csv'])
        f_freqs = open(f_freqs_path, 'r', newline='')
        print('Opening {0}'.format(f_freqs_path))
        f_velocities = open(f_velocities_path, 'r', newline='')
        print('Opening {0}'.format(f_velocities_path))
        freqs_reader = csv.reader(f_freqs, delimiter=',')
        velocities_reader = csv.reader(f_velocities, delimiter=',')
        freqs_data = []
        velocities_data = []

        for row in freqs_reader:
            freqs_data.append(row)
        for row in velocities_reader:
            velocities_data.append(row)

        if hasattr(self, 'numK') and hasattr(self, 'numBands'):
            numK = self.numK
            numBands = self.numBands
        else:
            self.numK = int(freqs_data[-1][0])
            self.numBands = len(freqs_data[0][5:])
            numK = self.numK
            numBands = self.numBands

        freqs_data.pop(0)  # pop out top row made up of list of strings only
        velocities_data.pop(0)
        f_freqs.close()
        f_velocities.close()

        # convert list of lists to appropriate numpy arrays. look at the list
        # of lists to understand the indexing
        k = np.zeros((numK, 3))
        kmag = np.zeros(numK)
        freqs = np.zeros((numK, numBands))
        vg = np.zeros((numK, 3, numBands))
        vgmag = np.zeros((numK, numBands))

        for i in range(0, numK):
            k[i, 0] = freqs_data[i][1]  # automatic typecast from str to float
            k[i, 1] = freqs_data[i][2]
            k[i, 2] = freqs_data[i][3]
            kmag[i] = freqs_data[i][4]
            freqs[i, :] = freqs_data[i][5:]
            # check to see is there is any actual data in velocities csv file
            if len(velocities_data) != 0:
                for j in range(0, numBands):
                    vg[i, :, j] = velocities_data[i][j+1].split()
                vgmag[i, :] = np.sqrt(vg[i, 0, :]**2 + vg[i, 1, :]**2 + vg[i, 2, :]**2)

        self.k = k
        self.kmag = kmag
        self.freqs = freqs
        if len(velocities_data) != 0:
            self.vg = vg
            self.vgmag = vgmag
        else:
            print('No velocity data in velocities csv file. Velocity attributes NOT created.')
        print('Parsing of csv files complete. New data attributes created!')


class h5Dataset:
    def __init__(self, mpb_h5_fobj, dset):
        self.dset = mpb_h5_fobj[dset]
        self.h5_fobj = mpb_h5_fobj

    def close(self):
        self.h5_fobj.close()


class EMField(h5Dataset):
    """
    Class for tensor components of the EM fields.
    Inherits methods of Class h5Dataset
    """
    def __init__(self, mpb_h5_fobj, mpbpostprocess=False):
        if not mpbpostprocess:
            self.xr = mpb_h5_fobj['x.r']
            self.xi = mpb_h5_fobj['x.i']
            self.yr = mpb_h5_fobj['y.r']
            self.yi = mpb_h5_fobj['y.i']
            self.zr = mpb_h5_fobj['z.r']
            self.zi = mpb_h5_fobj['z.i']
        else:
            self.xr = mpb_h5_fobj['x.r-new']
            self.xi = mpb_h5_fobj['x.i-new']
            self.yr = mpb_h5_fobj['y.r-new']
            self.yi = mpb_h5_fobj['y.i-new']
            self.zr = mpb_h5_fobj['z.r-new']
            self.zi = mpb_h5_fobj['z.i-new']
        self.h5_fobj = mpb_h5_fobj

    def create_complex(self):
        I = complex(0, 1)
        # [:] creates a copy of the dataset in RAM
        self.x = self.xr[:] + I*self.xi[:]
        self.y = self.yr[:] + I*self.yi[:]
        self.z = self.zr[:] + I*self.zi[:]
        print('Complex fields created. Closing .h5 file and deleting attributes')
        self.close()  # close .h5 files
        del self.xr
        del self.xi
        del self.yr
        del self.yi
        del self.zr
        del self.zi
        del self.h5_fobj


class epsilon_tensor(h5Dataset):
    """
    """
    def __init__(self, mpb_h5_fobj, mpbpostprocess=False):
        if not mpbpostprocess:
            self.xx = mpb_h5_fobj['epsilon.xx']
            self.xy = mpb_h5_fobj['epsilon.xy']
            self.yx = self.xy
            self.xz = mpb_h5_fobj['epsilon.xz']
            self.zx = self.xz
            self.yy = mpb_h5_fobj['epsilon.yy']
            self.yz = mpb_h5_fobj['epsilon.yz']
            self.zy = self.yz
            self.zz = mpb_h5_fobj['epsilon.zz']
        else:
            self.xx = mpb_h5_fobj['epsilon.xx-new']
            self.xy = mpb_h5_fobj['epsilon.xy-new']
            self.yx = self.xy
            self.xz = mpb_h5_fobj['epsilon.xz-new']
            self.zx = self.xz
            self.yy = mpb_h5_fobj['epsilon.yy-new']
            self.yz = mpb_h5_fobj['epsilon.yz-new']
            self.zy = self.yz
            self.zz = mpb_h5_fobj['epsilon.zz-new']

        self.h5_fobj = mpb_h5_fobj


def readfield(mpbobj=None, kindex=None, band=None, field_type=None,
              field_file=None, mpbpostprocess=False):
    """
    Reads the appropriate .h5 file from MPB output. Defaults to epsilon.h5

    Inputs
    ------
    mpbobj : Object of MPBBandStructure class
    kindex : an integer value used to read .h5 files
    band : an integer value used to read .h5 files
    field_type : 'e'
    field_file : read directly from .h5 file
    mpbpostprocess : False (defualt). mpb-data has not been run.
                    True. assumes the new datasets are embedded in the same .h5 file with -new append

    Returns
    -------

    Examples
    --------
    """

    # fname = '/home/nishan/Code/thales/MPB/w14/e.k03.b16.zeven.h5'
    if isinstance(mpbobj, MPBBandStructure):
        mpb = mpbobj
        path = mpb.path
        symm = mpb.symm
        if mpbobj.numK < 10:
            width = 1
        elif mpbobj.numK < 99:
            width = 2
    elif field_file:
        if not isinstance(field_file, str):
            print('field_file must be str')
            return None
    else:
        print('Input to readfield in invalid')
        return None


    if field_type == 'e' or field_type == 'h' or field_type == 'e.nonbloch.v' \
    or field_type == 'h.nonbloch.v':
        # open the .h5 file
        if field_file:
            f = h5py.File(field_file, 'r')
            field = EMField(f)
        else:
            if symm != '':
                fname = '/'.join([path, field_type + '.k' + str(kindex).zfill(width) + '.b' + str(band).zfill(width) + '.' + symm + '.h5'])
            else:
                 fname = '/'.join([path, field_type + '.k' + str(kindex).zfill(width) + '.b' + str(band).zfill(width) + '.h5'])
            f = h5py.File(fname, 'r')
            field = EMField(f, mpbpostprocess)
    elif field_type == 'epsilon':
        if field_file:
            f = h5py.File(field_file, 'r')
            field = epsilon_tensor(f, mpbpostprocess)
        else:
            fname = '/'.join([path, 'epsilon.h5'])
            f = h5py.File(fname, 'r')
            field = epsilon_tensor(f, mpbpostprocess)
    elif field_type == 'epsilon_isotropic_trace':
        if field_file:
            f = h5py.File(field_file, 'r')
            if not mpbpostprocess:
                print('Yo')
                field = h5Dataset(f, dset='data')
            else:
                field = h5Dataset(f, dset='data-new')
        else:
            fname = '/'.join([path, 'epsilon.h5'])
            f = h5py.File(fname, 'r')
            if not mpbpostprocess:
                field = h5Dataset(f, dset='data')
            else:
                field = h5Dataset(f, dset='data-new')
    else:
        print('MPBBandstructure:readfield:Invalid field_type entered')
        return None
    return field


def getscale(mpb, retstep=False):
    """
    Get x,y,z grids. ONLY FOR RECTANGULAR GRIDS DOES THIS MAKE SENSE
    """

    # read epsilon.h5 file for lattice vectors
    fname = '/'.join([mpb.path, 'epsilon.h5'])
    f = h5py.File(fname, 'r')
    epsilon = h5Dataset(f, dset='data')
    lattice_vecs_h5 = h5Dataset(f, dset='lattice vectors')
    lattice_vecs = lattice_vecs_h5.dset
    dim = len(epsilon.dset.shape)

    # Nx, Ny, Nz
    if dim == 1:
        Nx = epsilon.dset.shape[0]
        # For RECTANGULAR GRID ONLY
        (x,dx) = np.linspace(0, lattice_vecs[0,0], Nx, retstep=True)
    elif dim == 2:
        Nx = epsilon.dset.shape[0]
        Ny = epsilon.dset.shape[1]
        # For RECTANGULAR GRID ONLY
        (x,dx) = np.linspace(0, lattice_vecs[0,0], Nx, retstep=True)
        (y,dy) = np.linspace(0, lattice_vecs[1,1], Ny, retstep=True)
    elif dim == 3:
        Nx = epsilon.dset.shape[0]
        Ny = epsilon.dset.shape[1]
        Nz = epsilon.dset.shape[2]
        # For RECTANGULAR GRID ONLY
        (x,dx) = np.linspace(0, lattice_vecs[0,0], Nx, retstep=True)
        (y,dy) = np.linspace(0, lattice_vecs[1,1], Ny, retstep=True)
        (z,dz) = np.linspace(0, lattice_vecs[2,2], Nz, retstep=True)

    epsilon.close()

    if retstep:
        if dim == 1:
            return (x, dx)
        elif dim == 2:
            return (x, dx, y, dy)
        elif dim == 3:
            return (x, dx, y, dy, z, dz)
    else:
        if dim == 1:
            return (x)
        elif dim == 2:
            return (x, y)
        elif dim == 3:
            return (x, y, z)



# ks = []  # tuple of k points
# # look for line of the form 'solve_kpoint (0.3,0,0)'
# if re.match(r'solve_kpoint', line):
#     #  [(]\d[.]*[\d]*,\d[.]*[\d]*,\d[.]*[\d]*[)] selects any expressions of
#     # the form (0.1215,0.12,0.213) or (1,0,0.1564) etc.
#     #kpoint_match = re.search(r'[(]\d[.]*[\d]*,\d[.]*[\d]*,\d[.]*[\d]*[)]', line)
#     #k.append(kpoint_match.group())
#
#     # finds all values of the form 0.12315 and puts them in a list
#     (kx, ky, kz) = re.findall(r'\d[.]*[\d]*', line)
#     ks.append((float(kx), float(ky), float(kz)))  # each k-point is a tuple

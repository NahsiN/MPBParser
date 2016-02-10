# Parse the MPB output file for relevant quantities
# Inspired from Patterson's Matlab code

# parse the text file
import re
import csv
import numpy as np
import h5py

def strip_spaces(fields):
    ''' Strip spaces and newline characters from a list of strings '''
    # id(fields[:]) ≠ id(fields). fields[:] creates a new copy to allow
    # me to edit fields inside the loop (Remember Python passes by reference)
    # strip out spaces
    for element in fields[:]:
        fields.pop(0)
        fields.append(element.strip(' \n')) # \n is only needed for last element
    return fields

class MPBBandStructure:
    """
    Class for parsing MPB output
    """

    def __init__(self, out_file, symm):
        self.fname_mpb = out_file
        self.symm = symm
        #self.path = out_file.rstrip(re.search('[\w.]+.out', out_file).group())  # relies on file being named foo.out
        # Unix paths [0:-1] ensures not to include filename. Spaces allowed?
        self.path = '/'.join(out_file.split('/')[0:-1])

    def csvparser(self):
        """
        Parses MPB output file to two csv files while defining class data attributes
        """

        symm = self.symm
        fname_mpb = self.fname_mpb

        f_mpb = open(fname_mpb, 'r')
        f_freqs = open(symm+'freqs.csv', 'w', newline='')
        f_velocities = open(symm+'velocities.csv', 'w', newline='')
        freqs_writer = csv.writer(f_freqs)  # initialize a writer object
        velocities_writer = csv.writer(f_velocities)

        for line in f_mpb:
            # determine number of k points and bands
            if re.match(r'Computing [\d]* bands', line):
                self.numBands = int(re.findall(r'[\d]+', line)[0])
            if re.match(r'[\d]+ k-points', line):
                self.numK = int(re.findall(r'[\d]+', line)[0])

            # determine a header to write to the velocites .csv file
            # search for the line 'zevenfreqs:, k index, k1, .....'
            if re.search(r'k index', line):
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
        print('Two csv files {0} and {1} created.'.format(symm+'freqs.csv', symm+'velocities.csv'))
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
        numK = self.numK
        numBands = self.numBands

        f_freqs = open(symm+'freqs.csv', 'r', newline='')
        f_velocities = open(symm+'velocities.csv', 'r', newline='')
        freqs_reader = csv.reader(f_freqs, delimiter=',')
        velocities_reader = csv.reader(f_velocities, delimiter=',')
        freqs_data = []
        velocities_data = []

        for row in freqs_reader:
            freqs_data.append(row)
        for row in velocities_reader:
            velocities_data.append(row)
        freqs_data.pop(0)  # pop out top row made up of list of strings only
        velocities_data.pop(0)
        f_freqs.close()
        f_velocities.close()

        # convert list of lists to appropriate numpy arrays. look at the list of
        # lists to understand the indexing
        k = np.zeros((numK,3))
        kmag = np.zeros(numK)
        freqs = np.zeros((numK,numBands))
        vg = np.zeros((numK,3,numBands))
        vgmag = np.zeros((numK,numBands))

        for i in range(0,numK):
            k[i,0] = freqs_data[i][1]  # automatic typecast from str to float
            k[i,1] = freqs_data[i][2]
            k[i,2] = freqs_data[i][3]
            kmag[i] = freqs_data[i][4]
            freqs[i,:] = freqs_data[i][5:]
            for j in range(0,numBands):
                vg[i,:,j] = velocities_data[i][j+1].split()
            vgmag[i,:] = np.sqrt(vg[i,0,:]**2 + vg[i,1,:]**2 + vg[i,2,:]**2)

        self.k = k; self.kmag = kmag
        self.freqs = freqs
        self.vg = vg; self.vgmag = vgmag
        print('Parsing of csv files complete. New data attributes created!')

class EMField():
    def __init__(self, mpb_h5_fobj):
        self.xr = mpb_h5_fobj['x.r']
        self.xi = mpb_h5_fobj['x.i']
        self.yr = mpb_h5_fobj['y.r']
        self.yi = mpb_h5_fobj['y.i']
        self.zr = mpb_h5_fobj['z.r']
        self.zi = mpb_h5_fobj['z.i']
        self.h5_fobj = mpb_h5_fobj

    def close(self):
        self.h5_fobj.close()        

def readfield(kindex, band, mpb, field_type='e', path=''):
    #fname = '/home/nishan/Code/thales/MPB/w14/e.k03.b16.zeven.h5'
    if path == '':
        # choode path from mpb.root
        path = mpb.path
    if mpb.numK < 10:
        width = 1
    elif mpb.numK < 99:
        width = 2
    fname = '/'.join([path, field_type + '.k' + str(kindex).zfill(width) + '.b' + str(band).zfill(width) + '.' + mpb.symm + '.h5'])
    f = h5py.File(fname, 'r')
    field = EMField(f)
    return field

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

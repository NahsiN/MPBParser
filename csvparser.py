# Parse the MPB output file for relevant quantities
# Inspired from Patterson's Matlab code

# parse the text file
import re
import csv

# input parameteres
symm = 'zeven'
fname_mpb_output = '/home/nishan/Code/thales/MPB/w14/w14.out'

f_mpb = open(fname_mpb_output, 'r')
f_freqs = open(symm+'freqs.csv', 'w', newline='')
f_velocities = open(symm+'velocities.csv', 'w', newline='')
freqs_writer = csv.writer(f_freqs)  # initialize a writer object
velocities_writer = csv.writer(f_velocities)

for line in f_mpb:
    # write out freqs to .csv file
    # look for lines of the form 'zevenfreqs:,'
    if re.match(symm + r'freqs:,[\s\S]*', line):
        # create a list by splitting line
        fields = (line.lstrip(symm+'freqs:,')).split(',')
        # strip out spaces
        # id(fields[:]) ≠ id(fields). fields[:] creates a new copy to allow
        # me to edit fields inside the loop (Remember Python passes by reference)
        for element in fields[:]:
            fields.pop(0)
            fields.append(element.strip(' \n'))
        freqs_writer.writerow(fields)  # write out list as row in .csv file

    # write out velocities to .csv file
    if re.match(symm + r'velocity:,[\s\S]*', line):
        fields = (line.lstrip(symm+'velocity:,')).split(',')
        for element in fields[:]:
            fields.pop(0)
            fields.append(element.strip(' #()\n'))  # \n is only needed for last element
        velocities_writer.writerow(fields)  # write out list as row in .csv file
f_mpb.close()
f_freqs.close()
f_velocities.close()

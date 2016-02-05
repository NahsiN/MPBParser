# Parse the MPB output file for relevant quantities
# Inspired from Patterson's Matlab code

# parse the text file
import re
# get the number of k-points
f = open('/home/nishan/Code/thales/MPB/w14/w14.out', 'r')

ks = []  # tuple of k points
for line in f:
    # store the k-point
    if re.match(r'solve_kpoint', line):
        #  [(]\d[.]*[\d]*,\d[.]*[\d]*,\d[.]*[\d]*[)] selects any expressions of
        # the form (0.1215,0.12,0.213) or (1,0,0.1564) etc.
        #kpoint_match = re.search(r'[(]\d[.]*[\d]*,\d[.]*[\d]*,\d[.]*[\d]*[)]', line)
        #k.append(kpoint_match.group())

        # finds all values of the form 0.12315 and puts them in a list
        (kx, ky, kz) = re.findall(r'\d[.]*[\d]*', line)
        ks.append((float(kx), float(ky), float(kz)))  # each k-point is a tuple

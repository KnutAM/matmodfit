import sys
# Simple example of external simulation for matmodfit.
# Compares the two first material parameters to given values of E and nu
E = 211e3
nu= 0.33

# Check input options
if (sys.argv[1]=='1'):
    evec_out = True
else:
    evec_out = False

# Save file names
mpar_file = sys.argv[3]
error_file= sys.argv[4]
evec_file = sys.argv[5]

# Read in material parameters
mpar_fid = open(mpar_file, 'r')
mpar_str = mpar_fid.readlines()
mpar_fid.close()

mpar = [float(s.replace('\n', '')) for s in mpar_str]

# Calculate error vector and total error
evec = [(E-mpar[0])/E, (nu-mpar[1])/nu]
error = sum([e**2 for e in evec])

# Write error to file
error_fid = open(error_file, 'w')
error_fid.write('{:20.15e}'.format(error))
error_fid.close()

# Write error vector to file if requested
if evec_out:
    evec_fid = open(evec_file, 'w')
    for e in evec:
        evec_fid.write('{: 20.15e}'.format(e) + '\n')
    evec_fid.close()
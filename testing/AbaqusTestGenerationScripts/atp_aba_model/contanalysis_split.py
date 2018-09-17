import sys

def write_line(fid_out, line, epsz0, phi0, cstri0, cstro0, epsz_col, phi_col, cstri_col, cstro_col, h0, h1):
    the_line = line[:]
    # Convert strain values and write the line
    # Get actual strain values
    epsz  = float(the_line[epsz_col])
    phi   = float(the_line[phi_col])
    cstri = float(the_line[cstri_col])
    cstro = float(the_line[cstro_col])
    
    # Calculate strain values based on new reference lengths
    epsz1  = (epsz-epsz0)/(1+epsz0)
    phi1   = (phi-phi0)*h1/(h0*(1+epsz0))
    cstri1 = (cstri-cstri0)/(1+cstri0)
    cstro1 = (cstro-cstro0)/(1+cstro0)
    
    # Put new strain back into the_line
    the_line[epsz_col]  = '%16.8E'%epsz1
    the_line[phi_col]   = '%16.8E'%phi1
    the_line[cstri_col] = '%16.8E'%cstri1
    the_line[cstro_col] = '%16.8E'%cstro1
    
    # Write the split line to the new file
    fid_out.write(('%7s' + ' %16s'*(len(the_line)-1) + '\n')%tuple(the_line))
    

# Split the test data into several files, and modify the last to incorporate the changes simulating a new mounting of an extensometer
test_name = 'test_13'
h = [5.0, 4.0]     # Heights for each simulation, first much match the one in the Abaqus simulation
split_steps = [5]    # Which steps to split at the end of

# Column definitions (decrement with one due to pythons index from zero)
conv_index = [0, 3, 4, 6, 8, 10, 12, 14, 16, 18]    # Convert from full file to output file
# Strain indices in conv_index (note that python index from zero)
epsz_col  = 4-1
cstri_col = 8-1
cstro_col = 10-1
phi_col   = 6-1

exp_data  = test_name + '_raw.txt'         # Raw verification experimental data file (from matmodfit full sim)

# Get header:
fid_in      = open(exp_data, 'r')
the_string  = fid_in.readline()[1:]
the_title   = [the_string.split()[i] for i in conv_index]

# Set old strains to zero
epsz0  = 0.0
cstri0 = 0.0
cstro0 = 0.0
phi0   = 0.0
the_line = ['']

split_steps.append(float("inf"))

for i in range(len(h)):
    # Open new output file and write header/title line
    fid_out = open(test_name + '_exp_data_' + str(i+1) + '.txt', 'w')
    fid_out.write(('!%6s' + ' %16s'*(len(the_title)-1) + '\n')%tuple(the_title))
    
    if i>0: 
        # Need to write the initial starting line (i.e. the previous experiment last line)
        old_line[0]         = the_line[0]   # Change step to new step to avoid the zero step in the beginning
        write_line(fid_out, old_line, epsz0, phi0, cstri0, cstro0, 
                       epsz_col, phi_col, cstri_col, cstro_col, h[0], h[i])
        # Need also to write the first real line as this is already read
        write_line(fid_out, the_line, epsz0, phi0, cstri0, cstro0, 
                       epsz_col, phi_col, cstri_col, cstro_col, h[0], h[i])
    
    next_file = False
    while (not next_file):
        the_string = fid_in.readline()
        split_string = the_string.split()
        while(len(the_string)>0 and len(split_string)<len(conv_index)):
            the_string = fid_in.readline()
            split_string = the_string.split()
        try:
            old_line = the_line
            the_line = [split_string[i] for i in conv_index]
        except IndexError:
            if len(the_string)>0:   # Not due to end of file
                print('Something wrong? Line is:')
                print(the_string)
            # Either an error occurred or end of fid_in
            break
        step = int(the_line[0])
        if step>split_steps[i]:
            next_file = True
        else:
            write_line(fid_out, the_line, epsz0, phi0, cstri0, cstro0, 
                       epsz_col, phi_col, cstri_col, cstro_col, h[0], h[i])
    # End of while (not next_file)
    
    # Get the end values after split_step
    epsz0  = float(old_line[epsz_col])
    cstri0 = float(old_line[cstri_col])
    cstro0 = float(old_line[cstro_col])
    phi0   = float(old_line[phi_col])


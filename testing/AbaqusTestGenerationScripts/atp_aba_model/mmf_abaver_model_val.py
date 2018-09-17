import numpy as np

def setup_model(mtrl_model):
    if mtrl_model=='Chaboche_builtin' or mtrl_model=='Chaboche':
        #mpar= [   Emod,  nu,   sy0, Hiso, inv_kap_inf, Hk1, inv_b_inf1]
        mpar = [210.0e3, 0.3, 400.0, 10e3,      0.0025, 5e3,      0.002]
        nlgeom = False
        nstatv = 8
        umatpath = r'D:\BoxSync\knutan\PhD\CODING\Fortran_libraries\matmodfit\compiled_umats_abaqus\chaboche\chaboche.obj'
        umat_lib = 'chaboche'
    elif mtrl_model=='GFS':
        #mpar= [Gmod, Kmod,   Y0, Hiso, invYiso, delta, Hk1, invYk1, m1]
        mpar = [80e3, 175, 400.0, 10e3,  0.0025,   1.0, 5e3,  0.002, 1.0]
        nlgeom = True
        nstatv = 19
        umatpath = r'D:\BoxSync\knutan\PhD\CODING\Fortran_libraries\umats\GenFiniteStrain\Compiled\GFS_FE1_EV1_SC1\GFS_FE1_EV1_SC1-std.obj'
        umat_lib = 'gfs_af'
    else:
        print('Material model is not supported')
        return
    
    return mpar, nlgeom, nstatv, umatpath, umat_lib

def generate_amplitudes(total_time, dt):
    # Generate (triangular, sinusoidal, triangular and sinusoidal)
    # amplitude signals with frequencies (2, 4, 4, 2)/total_time
    time = np.linspace(0, total_time, int(total_time/dt))
    f1 = 2.0/total_time
    f2 = 4.0/total_time
    
    amp1 = triangular(f1*time)
    amp2 = np.sin(2*np.pi*f2*time)
    amp3 = triangular(f2*time)
    amp4 = np.sin(2*np.pi*f1*time)
    
    amp = []
    amp.append(list(zip(*[time,amp1])))
    amp.append(list(zip(*[time,amp2])))
    amp.append(list(zip(*[time,amp3])))
    amp.append(list(zip(*[time,amp4])))
    return amp
    
def triangular(time):
    return 2.0*(2.0*time - np.floor(2.0*time+0.5))*np.power(-1, np.floor(2.0*time+0.5))

# The input in this file can be change to adapt to different scenarios

# Input for matmodfit mesh type (only equidistant nodes allowed hence only nel, ri and ro specified instead of node_pos)
# Script automatically calculates equivalent Abaqus input
mmf_input_name = 'test_05'  # matmodfit input name (excluding suffix, '.inp' added automatically)
h0 = 5.0    
nel = 2     # Number of elements
ri = 0.0    # Inner radius
ro = 4.0    # Outer radius
node_pos = [ri+i*(ro-ri)/nel for i in range(nel+1)]
element_order = 2

# Always full integration during tests:
ngp = 2 if element_order==1 else 3
abaqus_bbar = (element_order==1)

# Choose material model, options are: 'Chaboche_builtin', 'Chaboche', 'GFS'
mtrl_model = 'GFS'
mpar, nlgeom, nstatv, umatpath, umat_lib = setup_model(mtrl_model)

# Loading (strain controlled)
total_time  = 2.00  # Total time (f0=2/T)
dt          = 0.002 # Fixed increment time
epsz_amp    = 0.04  # f = f0    (triangular)
gamma_amp   = 0.04  # f = 2*f0  (sinusoidal)
cstri_amp   = 0.01  # f = 2*f0  (triangular)
cstro_amp   = 0.02  # f = f0    (sinusoidal)
amp = generate_amplitudes(total_time, dt)
step_time = [total_time/8.0]*9
step_time[0] = step_time[0]/2.0
step_time[-1] = step_time[-1]/2.0

# matmodfit specific options (that doesn't affect Abaqus)
# Control type (Abaqus is run disp control (ctrl=1)), but should produce the same
# result also for different control modes. 
ctrl = [1, 1, 1, 1] # axial, torsion, inside, outside

# Generate file names
JOB_NAME = mmf_input_name + '_aba'
exp_data  = mmf_input_name + '_exp_data.txt'
mmf_input_file = mmf_input_name + '.inp'

# Convenience values
rm = (ri+ro)/2.0    # Center radially
hm = h0/2.0         # Center axially

# Constants
MTRL_NAME = 'MY_MATERIAL'
PART_NAME = 'MY_PART'
INST_NAME = 'MY_INSTANCE'
BOT_SET_NAME = 'BOT_SET'
TOP_SET_NAME = 'TOP_SET'
INS_SET_NAME = 'INS_SET'
OUT_SET_NAME = 'OUT_SET'

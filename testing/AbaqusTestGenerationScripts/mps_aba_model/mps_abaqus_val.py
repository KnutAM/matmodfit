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

    return 2.0*(2.0*time - np.floor(2.0*time+0.5))*np.power(-1, np.floor(2.0*time+0.5))

# Input for matmodfit mesh type (only equidistant nodes allowed hence only nel, ri and ro specified instead of node_pos)
# Script automatically calculates equivalent Abaqus input
mmf_input_name = 'test_16'  # matmodfit input name (excluding suffix, '.inp' added automatically)


# Choose material model, options are: 'Chaboche_builtin', 'Chaboche', 'GFS'
mtrl_model = 'GFS'
mpar, nlgeom, nstatv, umatpath, umat_lib = setup_model(mtrl_model)

# Loading (strain controlled)
total_time  = 2.00  # Total time (f0=2/T)
dt          = 0.005 # Fixed increment time
F11         = 1.5
F22         = 0.8
F33         = 1.05/(F11*F22)
F13         = 0.2
F23         = 1.5

# matmodfit specific options (that doesn't affect Abaqus)
# Control type (Abaqus is run disp control (ctrl=1)), but should produce the same
# result also for different control modes.

# Generate file names
JOB_NAME = mmf_input_name + '_aba'
exp_data  = mmf_input_name + '_exp_data.txt'
mmf_input_file = mmf_input_name + '.inp'


# Constants
MTRL_NAME = 'MY_MATERIAL'
PART_NAME = 'MY_PART'
INST_NAME = 'MY_INSTANCE'
BOT_SET_NAME = 'BOT_SET'
TOP_SET_NAME = 'TOP_SET'
XAX_SET_NAME = 'XAX_SET'
YAX_SET_NAME = 'YAX_SET'
ZAX_SET_NAME = 'ZAX_SET'

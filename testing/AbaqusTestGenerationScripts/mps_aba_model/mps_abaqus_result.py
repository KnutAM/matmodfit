from visualization  import *
from odbAccess      import *
import operator

import numpy as np

import mps_abaqus_val as val

def generate_mmf_input():
    write_exp_data()
    write_inputfile()
    
def write_exp_data():
    odb         = openOdb(path=(val.JOB_NAME+'.odb'))
    fid         = open(val.exp_data, 'w')
    
    if val.nlgeom:
        S = ['P11 [MPa]', 'P22 [MPa]', 'P33 [MPa]', 'P13 [MPa]', 'P23 [MPa]']
        E = ['F11-1.0', 'F22-1.0', 'F33-1.0', 'F13', 'F23']
    else:
        S = ['sig11 [MPa]', 'sig22 [MPa]', 'sig33 [MPa]', 'sig13 [MPa]', 'sig23 [MPa]']
        E = ['eps11', 'eps22', 'eps33', 'gam13', 'gam23']
    
    fid.write( ('!%6s' + ' %16s'*11 + '\n')% ('step', 'time [s]', 
                                             S[0], E[0], S[1], E[1],
                                             S[2], E[2], S[3], E[3],
                                             S[4], E[4]))
    
    tadd = 0.0
    stepnr = 0
    for step in odb.steps.values():
        stepnr = stepnr + 1
        tadd = tadd + write_step_output(fid, odb, step, stepnr, tadd)
    
    fid.close()
    
def write_step_output(fid, odb, step, stepnr, tadd):
    
    ## Deformation gradient values
    F11_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['F11_PT']
    F22_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['F22_PT']
    F33_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['F33_PT']
    F13_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['F13_PT']
    F23_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['F23_PT']
    
    F11 = get_average_displacement(step, F11_set, 'U1')
    F22 = get_average_displacement(step, F22_set, 'U2')
    F33 = get_average_displacement(step, F33_set, 'U3')
    F13 = get_average_displacement(step, F13_set, 'U1')
    F23 = get_average_displacement(step, F23_set, 'U2')
    
    ## Stress values
    P11_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['P11_PT']
    P22_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['P22_PT']
    P33_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['P33_PT']
    P13_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['P33_PT'] # TOP SURFACE, SAME AS P33
    P23_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets['P33_PT'] # TOP SURFACE, SAME AS P33
    
    P11 = get_reaction_force(step, P11_set, 'RF1')
    P22 = get_reaction_force(step, P22_set, 'RF2')
    P33 = get_reaction_force(step, P33_set, 'RF3')
    P13 = get_reaction_force(step, P13_set, 'RF1')
    P23 = get_reaction_force(step, P23_set, 'RF2')
    
    ## Time
    a_node = 'Node ' + F11_set.nodes[-1].instanceName + '.' + str(F11_set.nodes[-1].label)
    tmp    = step.historyRegions[a_node].historyOutputs['U1'].data
    time     = zip(*tmp)[0]
    
    ## Write output to file
    data = (time, P11, F11, P22, F22, P33, F33, P13, F13, P23, F23)
    data = zip(*data)
    # for t in time:
        # fid.write(('%7.0F \n') % t)
    
    
    for runt, P11, F11, P22, F22, P33, F33, P13, F13, P23, F23 in data:
        fid.write(('%7.0F' + ' %16.5E'*11 + '\n') % (stepnr, runt+tadd, 
                                                     P11, F11, P22, F22, 
                                                     P33, F33, 
                                                     P13, F13, P23, F23))
    return time[-1]
     
def get_average_displacement(the_step, the_node_set, the_field):
#   multiadd = lambda a,b: map(operator.add, a, b)
    Nnod = len(the_node_set.nodes)
    for i in range(Nnod):
        
        the_node = 'Node ' + the_node_set.nodes[i].instanceName + '.' + str(the_node_set.nodes[i].label)
        the_data = the_step.historyRegions[the_node].historyOutputs[the_field].data
        
        temp = zip(*the_data)[1]                # Extract the reaction values
        if i==0:
            out_data = [0.0]*len(temp)          
            
        out_data = multimean(out_data, temp, Nnod)     # Add contribution from current node to each time increment
    
    
    return out_data
     
def get_reaction_force(the_step, the_node_set, the_field):
#   multiadd = lambda a,b: map(operator.add, a, b)
    
    for i in range(len(the_node_set.nodes)):
        
        the_node = 'Node ' + the_node_set.nodes[i].instanceName + '.' + str(the_node_set.nodes[i].label)
        the_data = the_step.historyRegions[the_node].historyOutputs[the_field].data
        
        temp = zip(*the_data)[1]                # Extract the reaction values
        if i==0:
            out_data = [0.0]*len(temp)          
            
        out_data = multiadd(out_data, temp)     # Add contribution from current node to each time increment
    
    return out_data
    
def multiadd(a,b):
    N = len(a)
    c = [0.0]*N
    for i in range(N):
        c[i] = a[i]+b[i]
    
    return c
    
def multimean(a,b,NUM):
    N = len(a)
    c = [0.0]*N
    for i in range(N):
        c[i] = a[i]+b[i]/NUM
    
    return c
    
    
def write_inputfile():
    fid = open(val.mmf_input_file, 'w')
    ## Write some header information
    fid.write('! This input file was automatically generated by write_inputfile() in mps_abaqus_result.py\n')
    fid.write('! The input file is used for verifying that matmodfit\'s mps produce the same output as Abaqus\n')
    fid.write('! Ramp loading with F11 = ' + str(val.F11) + ', F22 = ' + str(val.F22) + ', F33 = ' + str(val.F33) + 
              ', F13 = ' + str(val.F13) + ', F23 = ' + str(val.F23) + ', F12=F21=F31=F32=0.0\n')
    
    if val.nlgeom:
        #[11,22,33,12,23,31,13,21,32]
        exp_info = '1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 11, 12, 0, 0, 9, 10, 0, 0, 0, 0, 0'
        ctrl     = list_2_str([1,1,1,1,1,1,1,1,1])
    else:
        #[11,22,33,12,13,23]
        exp_info = '1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 9, 10, 11, 12, 0'
        ctrl     = list_2_str([1,1,1,1,1,1])
    
    ## Write global settings
    fid.write(' *run_type\n' + '1\n')
    fid.write(' *umat_lib\n' + val.umat_lib + '\n')
    fid.write(' *param_init\n' + list_2_str(val.mpar) + '\n')
    fid.write(' *nstatv\n' + str(val.nstatv) + '\n')
    fid.write(' *nlgeom\n' + bool_2_intstr(val.nlgeom) + '\n')
    fid.write('\n') # Extra spacing between categories to facilitate reading
    
    ## Write simulation settings
    fid.write('<Simulation>\n' + '2\n')
    fid.write('<<outp>>\n')
    fid.write(' *result_inclexp\n' + '1\n')
    # Error settings
    fid.write('<<err>>\n')
    fid.write(' *err_norm_met\n' + '2\n')
    # Load and exp_data
    fid.write('<<exp>>\n')
    fid.write(' *ctrl\n' + '1\n' + '1, ' + ctrl + '\n')
    fid.write(' *exp_data\n' + val.exp_data + '\n')
    
    fid.write(' *exp_info\n' + exp_info + '\n')
    # Iteration settings
    fid.write('<<iter>>\n')
    fid.write(' *time_incr\n' + '1\n' + '1, ' + str(val.dt))
    
    fid.close()
    
def list_2_str(the_list):
    out_str = str(the_list)
    return out_str[1:-1]    # Remove brackets
    
def bool_2_intstr(the_bool):
    return '1' if the_bool else '0'
    
if __name__ == '__main__':              
    generate_mmf_input()                              # run the main function
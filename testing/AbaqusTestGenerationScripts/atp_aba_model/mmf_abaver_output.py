from visualization  import *
from odbAccess      import *
import operator

import numpy as np

import mmf_abaver_model_val as val

def generate_mmf_input():
    write_exp_data()
    write_inputfile()
    
def write_exp_data():
    odb         = openOdb(path=(val.JOB_NAME+'.odb'))
    fid         = open(val.exp_data, 'w')
    def_geom_fid= open('def_geom.txt', 'w')
    fid.write( ('!%6s' + ' %16s'*9 + '\n') % ('step', 'time [s]', 
                                             'axial force [N]', 'eps_z [-]',
                                             'torque [Nmm]', 'rotation [rad]', 
                                             'p_i [MPa]', 'cstr_i [-]',
                                             'p_o [MPa]', 'cstr_o [-]'))
    
    tadd = 0.0
    stepnr = 0
    for step in odb.steps.values():
        stepnr = stepnr + 1
        tadd = tadd + write_step_output(fid, odb, step, stepnr, tadd, def_geom_fid)
    
    fid.close()
    def_geom_fid.close()
    
def write_step_output(fid, odb, step, stepnr, tadd, def_geom_fid):
    ## Extract node sets and nodes
    top_node_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets[val.TOP_SET_NAME]
    ins_node_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets[val.INS_SET_NAME]
    out_node_set  = odb.rootAssembly.instances[val.INST_NAME].nodeSets[val.OUT_SET_NAME]
    
    top_node = 'Node ' + top_node_set.nodes[-1].instanceName + '.' + str(top_node_set.nodes[-1].label)
    ins_node = 'Node ' + ins_node_set.nodes[-1].instanceName + '.' + str(ins_node_set.nodes[-1].label)
    out_node = 'Node ' + out_node_set.nodes[-1].instanceName + '.' + str(out_node_set.nodes[-1].label)
    
    ## Get displacement values
    rota     = step.historyRegions[top_node].historyOutputs['UR2'].data
    delta_z  = step.historyRegions[top_node].historyOutputs['U2'].data
    delta_ri = step.historyRegions[ins_node].historyOutputs['U1'].data
    delta_ro = step.historyRegions[out_node].historyOutputs['U1'].data
    
    rota     = zip(*rota)   #"Transpose" tuple
    time     = rota[0]
    rota     = rota[1]
    delta_z  = zip(*delta_z)[1]
    delta_ri = zip(*delta_ri)[1]
    delta_ro = zip(*delta_ro)[1]
    
    # Convert to strains
    eps_z = [dz/val.h0 for dz in delta_z]
    cstr_i = [ (dri/val.ri if (val.ri>1e-10) else 0.0) for dri in delta_ri]
    cstr_o = [dro/val.ro for dro in delta_ro]
    
    ## Get reaction forces
    torque  = get_reaction_force(step, top_node_set, 'RM2')
    aforce  = get_reaction_force(step, top_node_set, 'RF2')
    iforce  = get_reaction_force(step, ins_node_set, 'RF1')
    oforce  = get_reaction_force(step, out_node_set, 'RF1')
    
    # Calculate pressures
    p_i = []
    p_o = []
    for dz, dri, dro, f_i, f_o in zip(*(delta_z, delta_ri, delta_ro, iforce, oforce)):
        if val.nlgeom:
            h = val.h0 + dz
            ri = val.ri + dri
            ro = val.ro + dro
        else:
            h = val.h0
            ri = val.ri
            ro = val.ro
            
        p_i.append(+f_i/(2.0*np.pi*ri*h) if (val.ri>1e-10) else 0.0)
        p_o.append(-f_o/(2.0*np.pi*ro*h))
    
    ## Write output to file
    data = (time, aforce, eps_z, torque, rota, p_i, cstr_i, p_o, cstr_o)
    data = zip(*data)
    
    for runt, forc, epsz, torq, rota, p_i, cstr_i, p_o, cstr_o in data:
        fid.write(('%7.0F' + ' %16.5E'*9 + '\n') % (stepnr, runt+tadd, 
                                                     forc, epsz, 
                                                     torq, rota, 
                                                     p_i, cstr_i, 
                                                     p_o, cstr_o))
    
    # Write out to file with updated node_pos
    def_node_pos = []
    for node in top_node_set.nodes:
        node_name = 'Node ' + node.instanceName + '.' + str(node.label)
        node_data = step.historyRegions[node_name].historyOutputs['U1'].data
        def_node_pos.append(node.coordinates[0] + node_data[-1][1])
    
    def_node_pos.sort()
    def_geom_fid.write('Step ' + str(stepnr) + ': deformed node_pos: ')
    def_geom_fid.write(('%16.5E'*len(def_node_pos) + '\n')%tuple(def_node_pos))
    
    return time[-1]
     
def get_reaction_force(the_step, the_node_set, the_field):
    multiadd = lambda a,b: map(operator.add, a, b)
    
    for i in range(len(the_node_set.nodes)):
        
        the_node = 'Node ' + the_node_set.nodes[i].instanceName + '.' + str(the_node_set.nodes[i].label)
        the_data = the_step.historyRegions[the_node].historyOutputs[the_field].data
        
        temp = zip(*the_data)[1]                # Extract the reaction values
        if i==0:
            out_data = [0.0]*len(temp)          
            
        out_data = multiadd(out_data, temp)     # Add contribution from current node to each time increment
    
    return out_data
    
    
def write_inputfile():
    fid = open(val.mmf_input_file, 'w')
    ## Write some header information
    fid.write('! This input file was automatically generated by write_inputfile() in mmf_abaver_output.py\n')
    fid.write('! The input file is used for verifying that matmodfit produce the same output as Abaqus\n')
    fid.write(('! amplitudes: (epsz, phi, cstri, cstro): (' + ' %0.6F '*4 + ')\n')%(val.epsz_amp, val.gamma_amp, val.cstri_amp, val.cstro_amp))
    
    ## Write global settings
    fid.write(' *run_type\n' + '1\n')
    fid.write(' *umat_lib\n' + val.umat_lib + '\n')
    fid.write(' *param_init\n' + list_2_str(val.mpar) + '\n')
    fid.write(' *nstatv\n' + str(val.nstatv) + '\n')
    fid.write(' *nlgeom\n' + bool_2_intstr(val.nlgeom) + '\n')
    fid.write('\n') # Extra spacing between categories to facilitate reading
    
    ## Write simulation settings
    fid.write('<Simulation>\n' + '1\n')
    fid.write('<<outp>>\n')
    fid.write(' *result_inclexp\n' + '1\n')
    # Mesh and geometry
    fid.write('<<mesh1d>>\n')
    fid.write(' *h0\n' + str(val.h0) + '\n')
    fid.write(' *node_pos\n' + list_2_str(val.node_pos) + '\n')
    fid.write(' *ngp\n' + str(val.ngp) + '\n')
    fid.write(' *element_order\n' + str(val.element_order) + '\n')
    fid.write(' *abaqus_bbar\n' + bool_2_intstr(val.abaqus_bbar) + '\n')
    # Error settings
    fid.write('<<err>>\n')
    fid.write(' *err_norm_met\n' + '2\n')
    # Load and exp_data
    fid.write('<<exp>>\n')
    fid.write(' *ctrl\n' + '1\n' + '1, ' + list_2_str(val.ctrl) + '\n')
    fid.write(' *exp_data\n' + val.exp_data + '\n')
    fid.write(' *exp_info\n' + '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0\n')
    # Iteration settings
    fid.write('<<iter>>\n')
    fid.write(' *time_incr\n' + '1\n' + '1, ' + str(val.dt))
    
    fid.close()
    
def list_2_str(the_list):
    out_str = str(the_list)
    return out_str[1:-1]    # Remove brackets
    
def bool_2_intstr(the_bool):
    return '1' if the_bool else '0'
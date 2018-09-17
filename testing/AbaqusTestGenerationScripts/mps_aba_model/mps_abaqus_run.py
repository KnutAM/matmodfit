import time
from abaqusConstants import *
from mps_abaqus_setup import setup_model
from mps_abaqus_result import generate_mmf_input

start_time = time.time()
fid = open('py_aba_log.log', 'w')

verif_job = setup_model()
setup_time = time.time() - start_time
fid.write("Setup time = %0.2F s\n" % setup_time)

verif_job.submit()
verif_job.waitForCompletion()
simulation_time = time.time() - start_time - setup_time
fid.write("Simulation time = %0.2F s\n" % simulation_time)

generate_mmf_input()
output_time = time.time() - start_time - setup_time - simulation_time
fid.write("Output time = %0.2F s\n" % output_time)

fid.write("Total time = %0.2F s\n" % (time.time() - start_time))
fid.close()
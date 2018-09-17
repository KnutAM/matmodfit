from __future__ import print_function
import sys
import os
import shutil
import platform
        

def main(input_args):
    # Input argument is output directory for models
    if len(input_args)>1:
        output_dir = input_args[1]
        output_dir = output_dir + ('' if output_dir[-1]=='/' else '/')
    else:
        output_dir = '../build_examples/'
        
    if len(input_args)>2:
        linux_compiler = input_args[2]
    else:
        linux_compiler = 'ifort'
        
    

    if platform.system()=='Linux':
        dll_suffix  = '.so'
        dll_location= '' 
        cmake_setup = 'FC=' + linux_compiler + ' cmake ../src'
    elif platform.system()=='Windows':
        dll_suffix = '.dll'
        dll_location = 'Release/'
        cmake_setup = 'cmake -A "x64" ../src'
    elif platform.system()=='Darwin':
        print('MAC not supported by this script')
    else:
        print('Unknown platform (os), not supported by this script')
    
    # Create output directory if needed
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    # Build user simulation example
    example_dir = '../examples/user_simulation/'
    output_name = 'usim_example'
    bld_example(cmake_setup, dll_suffix, dll_location, example_dir, output_dir, output_name)
    
    # Build user optimization example
    example_dir = '../examples/user_optimization/'
    output_name = 'uopt_example'
    bld_example(cmake_setup, dll_suffix, dll_location, example_dir, output_dir, output_name)
    
    # Build user material parameter scaling example
    example_dir = '../examples/user_mpar_scaling/'
    output_name = 'user_mpar_scale_example'
    bld_example(cmake_setup, dll_suffix, dll_location, example_dir, output_dir, output_name)
    
    

def bld_example(cmake_setup, dll_suffix, dll_location, example_dir, output_dir, output_name):
    # cmake_setup:  String with the cmake command to set up the build environment
    #               E.g. cmake ../src
    # dll_suffix:   '.dll' on windows and '.so' on linux
    # example_dir:  Rel or abs path to the directory containing the src folder for the example
    # output_dir:   Rel or abs path to the directory to which the built example library should be copied to
    # output_name:  Name of the resulting file
    
    old_dir = os.getcwd() + '/'
    os.chdir(example_dir)
    
    
    build_dir = 'build/'
    
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)
    
    os.mkdir(build_dir)
    os.chdir(build_dir)
    
    os.system(cmake_setup)
    os.system('cmake --build . --config RELEASE')
    
    os.chdir(old_dir)   # Change to old directory to allow output_dir to have a relative path
    shutil.copy(example_dir + build_dir + dll_location + output_name + dll_suffix, output_dir)
    os.chdir(example_dir + build_dir)

    # Change back to old directory before exiting function
    os.chdir(old_dir)
    

if __name__ == '__main__':              
    main(sys.argv)                              # run the main function

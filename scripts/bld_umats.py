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
        output_dir = '../build_umats/'
        
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
    
    # Build general small strain models:
    model_dir = '../umat/GeneralSmallStrainPlasticity/'
    bld_gss(cmake_setup, dll_suffix, dll_location, model_dir, output_dir)
    
    # Build Meyer2018a models
    model_dir = '../umat/Meyer2018a/'
    bld_meyer2018a(cmake_setup, dll_suffix, dll_location, model_dir, output_dir)
    

def bld_gss(cmake_setup, dll_suffix, dll_location, model_dir, output_dir):
    # cmake_setup:  String with the cmake command to set up the build environment
    #               E.g. cmake ../src
    # dll_suffix:   '.dll' on windows and '.so' on linux
    # model_dir:    Rel or abs path to the directory containing the src folder for the model
    # output_dir:   Rel or abs path to the directory in which the built models should be copied to
    
    old_dir = os.getcwd() + '/'
    os.chdir(model_dir)
    
    base_models = ['chaboche', 'delobelle', 'ohnowang']
    rate_models = ['norate', 'norton', 'cowsym', 'delobelle']
    
    build_dir = 'build/'
    
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)
    
    os.mkdir(build_dir)
    os.chdir(build_dir)
    
    os.system(cmake_setup)
    
    for bm in base_models:
        for rm in rate_models:
            modstr = ' -D basemod=' + bm + ' -D ratemod=' + rm
            os.system('cmake' + modstr + ' .')
            os.system('cmake --build . --config RELEASE')
            if rm=='norate':
                modname = bm + dll_suffix
            else:
                modname = bm + '_' + rm + dll_suffix
                
            os.chdir(old_dir)   # Change to old directory to allow output_dir to have a relative path
            shutil.copy(model_dir + build_dir + dll_location + modname, output_dir)
            os.chdir(model_dir + build_dir)

    # Change back to old directory before exiting function
    os.chdir(old_dir)
    
def bld_meyer2018a(cmake_setup, dll_suffix, dll_location, model_dir, output_dir):
    # cmake_setup:  String with the cmake command to set up the build environment
    #               E.g. cmake ../src
    # dll_suffix:   '.dll' on windows and '.so' on linux
    # model_dir:    Rel or abs path to the directory containing the src folder for the model
    # output_dir:   Rel or abs path to the directory in which the built models should be copied to
    
    old_dir = os.getcwd() + '/'
    os.chdir(model_dir)
    
    models = ['gfs_bc', 'gfs_ow']
    
    build_dir = 'build/'
    
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)
    
    os.mkdir(build_dir)
    os.chdir(build_dir)
    
    os.system(cmake_setup)
    
    for mod in models:
        modstr = ' -D model=' + mod
        os.system('cmake' + modstr + ' .')
        os.system('cmake --build . --config RELEASE')
        modname = mod + dll_suffix
        
        os.chdir(old_dir)   # Change to old directory to allow output_dir to have a relative path
        shutil.copy(model_dir + build_dir + dll_location + modname, output_dir)
        os.chdir(model_dir + build_dir)
    
    # Change back to old directory before exiting function
    os.chdir(old_dir)

if __name__ == '__main__':              
    main(sys.argv)                              # run the main function

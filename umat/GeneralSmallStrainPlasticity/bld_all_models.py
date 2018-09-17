import sys
import os
import shutil




def main(input_args):
    # Input is cmake additional arguments
    # E.g. build on windows with Visual studio 2012 64 bit:
    # python build_all.py "-G \"Visual Studio 11 2012 Win64\""
    base_models = ['chaboche', 'delobelle', 'ohnowang']
    rate_models = ['norate', 'norton', 'cowsym', 'delobelle']
    
    builddir = 'build/'
        
    if len(input_args)>1:
        cmake_addarg = input_args[1]
    else:
        cmake_addarg = ''
    
    if os.path.exists(builddir):
        shutil.rmtree(builddir)
    
    os.mkdir(builddir)
    os.chdir(builddir)
    
    os.system('cmake ' + cmake_addarg + ' ../src')
    
    for bm in base_models:
        for rm in rate_models:
            modstr = ' -D basemod=' + bm + ' -D ratemod=' + rm
            os.system('cmake' + modstr + ' .')
            os.system('cmake --build . --config RELEASE')

if __name__ == '__main__':              
    main(sys.argv)                              # run the main function

import sys
import time
import os
import subprocess
import shutil

def run_cases(input_args):
    tmpdir   = 'tmp_delete_before_commit'
    srcdir   = 'TestCaseFolder'
    resfile  = 'test_report.txt'
    casefile = 'test_cases.txt'
    
    # Create tmp folder to run test cases in (delete folder if exists)
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
        
    shutil.copytree(srcdir, tmpdir)
    shutil.copy(casefile, tmpdir)
    
    
    os.chdir(tmpdir)
    

    cases = get_cases(casefile)
    
    result_fid = open(resfile, 'w')
    
    if len(input_args)>1:
        try: 
            casenr = [int(a) for a in input_args[1:]]
        except ValueError:
            print('Only numeric input of the case rows in ' + casefile + ' is allowed')
            return
        if any([c<=0 for c in casenr]) or any([c>len(cases) for c in casenr]):
            print('Row nr ' + str(casenr) + ' is not defined in ' + casefile)
            return
        else:
            cases = [cases[c-1] for c in casenr]
    
    for c in cases:
        res = run_mmf(c['inp'], c['max_time'])
        write_result(c, res, result_fid)
        
    result_fid.close()
    shutil.copy(resfile, '..')
    os.chdir('..')


def run_mmf(simname, maxtime):
    result = {}
    
    start_time = time.time()
    try:
        sp_result = subprocess.run(['matmodfit', simname + '.inp'], 
                                   timeout=maxtime,
                                   stdout = subprocess.PIPE)
        result['output']  = sp_result.stdout.decode("utf-8") 
        result['success'] = sp_result.returncode == 0
        result['timeout'] = False
    except subprocess.TimeoutExpired as timeout:
        result['timeout'] = True
        result['output'] = timeout.output.decode("utf-8") 
        result['success'] = False

    result['time'] = time.time() - start_time
    result['error'] = 0.0
    
    if result['success']:
        try:
            fid = open(simname + '.res', 'r')
            fid.readline()    # Skip headers
            fid.readline()    # Skip headers
            result['error'] = float(fid.readline())    # 
        except FileNotFoundError:
            result['success'] = False
        except ValueError:
            result['success'] = False
        else:
            fid.close()

    return result
    
    
def get_cases(case_file):
    fid = open(case_file, 'r')
    fid.readline()  # Read first header line (scrap)
    tmp = fid.readline()  # Read case line
    data = tmp.split()
    cases = []
    while len(data)>=4:
        cases.append({'id': int(data[0]), 'inp': data[1], 'max_error': float(data[2]), 'max_time': float(data[3])})
        tmp = fid.readline()
        data = tmp.split()
    
    fid.close()
    return cases

    
def write_result(case, result, fid):
    fid.write('Case ' + str(case['id']) + ': ')
    if result['success'] and result['error']<case['max_error']:
        print('Case ' + str(case['id']) + ': PASS')
        fid.write('PASS\n')
        fid.write('   error = {:.2E}\n'.format(result['error']))
        fid.write('   time  = {:.2E} s\n'.format(result['time']))
    elif result['success']:
        print('Case ' + str(case['id']) + ': FAIL (ERROR = ' + str(result['error']) + ')')
        fid.write('FAIL (error too large)\n')
        fid.write('   error = {:.2E}\n'.format(result['error']))
        fid.write('   time  = {:.2E} s\n'.format(result['time']))
    elif result['timeout']:
        print('Case ' + str(case['id']) + ': FAIL (too much time spent)')
        fid.write('FAIL (timeout)\n')
        fid.write('simulation output: \n')
        fid.write(result['output'])
    else:
        print('Case ' + str(case['id']) + ': FAIL (program did not complete successfully)')
        fid.write('FAIL (unsuccessful program execution)\n')
        fid.write('simulation output: \n')
        fid.write(result['output'])
    
if __name__ == '__main__':              
    run_cases(sys.argv)                              # run the main function

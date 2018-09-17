import os
import shutil


tmpdir   = 'tmp_delete_before_commit'
resfile  = 'test_report.txt'


if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
        
try:
    os.remove(resfile)
except OSError:
    pass
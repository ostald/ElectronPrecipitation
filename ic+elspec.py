import os
import ic4elspec
import time
import sys

matlabroot_dir = "/Applications/MATLAB_R2022b.app/bin/./matlab"
if sys.platform == 'linux2':
    matlabroot_dir = "/usr/local/bin/matlab"

cwd = os.getcwd()
print('Current working Dir: ', cwd)

call_matlab = matlabroot_dir + " -sd \"" + cwd + "/ELSPEC-2022\" -batch "

#check if plotting folder exists:
if not os.path.isdir('log/testing/plots'):
    print('Creating log/testing/plots')
    os.mkdir('log/testing/plots')

#start model
i = 0
start_t = time.time()
os.system(call_matlab + "\"ElSpec_IC\" -nodisplay")
ic4elspec.ic('log/testing/', 'ElSpec-iqt_IC_', i)

i = 1
print('First Iteration:', time.time() - start_t, 's')

while True:
    os.system(call_matlab + "\"ElSpec_IC_iter("+str(i)+")\" -nodisplay")
    ic4elspec.ic('log/testing/',
                 'ElSpec-iqt_IC_', i)
    i = i+1
    print('Mean Iteration Duration:', (time.time() - start_t)/i, 's')

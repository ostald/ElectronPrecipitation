import os
import ic4elspec
import time

matlabroot_dir = "/Applications/MATLAB_R2022b.app/bin/./matlab"
cwd = os.getcwd()
print('Current working Dir: ', cwd)

#check if plotting folder exists:
if not os.path.isdir('log/testing/plots'):
    print('Creating log/testing/plots')
    os.mkdir('log/testing/plots')

#start model
i = 0
#os.system(matlabroot_dir + " -sd \"" + cwd + "/ELSPEC-2022\" "+
#          "-batch \"ElSpec_IC\" -nodisplay")

#ic4elspec.ic('log/testing/',
#             'ElSpec-iqt_IC_', i)
i = 1

while True:
    osstr = matlabroot_dir + " -sd \"" + cwd + "/ELSPEC-2022\" "+\
            "-batch \"ElSpec_IC_iter("+str(i)+")\" -nodisplay"
    os.system(osstr)
    ic4elspec.ic('log/testing/',
                 'ElSpec-iqt_IC_', i)
    i = i+1 




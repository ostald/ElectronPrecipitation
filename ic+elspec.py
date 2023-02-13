import os
import ic4elspec
import time
import sys
import setup

matlabroot_dir = "/Applications/MATLAB_R2022b.app/bin/./matlab"
if sys.platform == 'linux2':
    matlabroot_dir = "/usr/local/bin/matlab"

cwd = os.getcwd()
print('Current working Dir: ', cwd)


#defining file paths
msis_config = cwd + '/Data/other/msis.rtf'
iri_config  = cwd + '/Data/other/iri.txt'
chemistry_config = cwd + '/Data/other/Reaction rates full set.txt'
path_eiscat_data = cwd + '/Data/Eiscat'

mixf = 0

#setup log directory
setup_ = setup.setup(msis_config, iri_config, chemistry_config, path_eiscat_data, mixf = mixf)#, no_timestamp=True)
log_directory = setup_._log_directory

call_matlab = matlabroot_dir + " -sd \"" + cwd + "/ELSPEC-2022\" -batch "

#start model
i = 0
start_t = time.time()
os.system(call_matlab + "\"ElSpec_IC(\\\""+log_directory+"\\\")\" -nodisplay")

while True:
    ic4elspec.ic(log_directory, 'ElSpec-iqt_IC_', i, mixf = mixf)
    i = i+1
    print('Mean Iteration Duration:', (time.time() - start_t)/i, 's')
    os.system(call_matlab + "\"ElSpec_IC_iter("+str(i)+", \\\""+log_directory+"\\\")\" -nodisplay")

import os
import ic4elspec
import time
import sys
import setup
import loadmat


#defining file paths
msis_config = '/Users/ost051/Documents/PhD/Electron Precipitation/Data/other/msis.rtf'
iri_config  = '/Users/ost051/Documents/PhD/Electron Precipitation/Data/other/iri.txt'
chemistry_config = '/Users/ost051/Documents/PhD/Electron Precipitation/Data/other/Reaction rates full set.txt'
path_eiscat_data = '/Users/ost051/Documents/PhD/Electron Precipitation/Data/Eiscat'

mixf = 1

#setup log directory
setup_ = setup.setup(msis_config, iri_config, chemistry_config, path_eiscat_data, mixf)#, no_timecode=True)
log_directory = setup_._log_directory

matlabroot_dir = "/Applications/MATLAB_R2022b.app/bin/./matlab"
if sys.platform == 'linux2':
    matlabroot_dir = "/usr/local/bin/matlab"

cwd = os.getcwd()
print('Current working Dir: ', cwd)

#call_matlab = matlabroot_dir + " -sd \"" + cwd + "/ELSPEC-2022\" -r "
call_matlab = matlabroot_dir + " -sd \"" + cwd + "/ELSPEC-2022\" -batch "

#start model
i = 0
start_t = time.time()
os.system(call_matlab + "\"ElSpec_IC(\\\""+log_directory+"\\\")\" -nodisplay")

while True:
    ic4elspec.ic(log_directory, chemistry_config, 'ElSpec-iqt_IC_', i, mixf = mixf)
    i = i+1
    print('Mean Iteration Duration:', (time.time() - start_t)/i, 's')
    os.system(call_matlab + "\"ElSpec_IC_iter("+str(i)+", \\\""+log_directory+"\\\")\" -nodisplay")

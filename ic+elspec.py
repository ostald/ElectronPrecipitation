import os
import ic4elspec
import time
import sys
import setup
import loadmat

dirname = 'test_SlowElSpec'

#defining file paths
msis_config = '../Data/other/msis.rtf'
iri_config  = '../Data/other/iri.txt'
chemistry_config = '../Data/other/Reaction rates full set.txt'
#path_eiscat_data = '/mnt/data/bjorn/EISCAT/Analysed/2012-12-11_arc1_4@uhf' #Data/Eiscat'
#path_eiscat_data = '../Data/2012-12-11_arc1_4@uhf' #Data/Eiscat'
#path_eiscat_data = '../Data/Eiscat/pp'
path_eiscat_data = '../Data/Eiscat/pp'
result_path = "../Results"

mixf = 0

#setup log directory
setup_ = setup.setup(result_path, msis_config, iri_config, chemistry_config, path_eiscat_data, mixf, dirname=dirname)#, no_timecode=True)
log_directory = setup_._log_directory
print(log_directory)

matlabroot_dir = "/Applications/MATLAB_R2023b.app/bin/./matlab"
if sys.platform == 'linux':
    matlabroot_dir = "/usr/local/bin/matlab"

cwd = os.getcwd()

#call_matlab = matlabroot_dir + " -sd \"" + cwd + "/ELSPEC-2022\" -r "
#call_matlab = matlabroot_dir + " -sd \"" + cwd + "/ELSPEC-2022\" -batch "
call_matlab = matlabroot_dir + " -sd \"" + cwd + "/../ELSPEC\" -batch "

#start model
i = 0
start_t = time.time()
#os.system(call_matlab + "\"ElSpec_IC(\\\""+log_directory+"\\\")\" -nodisplay")

# import datetime
# btime = [2006, 12, 12, 19, 30, 00]
# etime = [2006, 12, 12, 22, 00, 00]
#
# bdt = datetime.datetime(*btime)
# edt = datetime.datetime(*etime)
#
# stime = [bdt]
# while stime[-1] < edt:
#     stime = [*stime, stime[-1] + datetime.timedelta(minutes=15)]
#
# for it in range(1, len(stime)-1):
#     if it == 0:
#         #no ne_init
#     #call matlab(btime = stime[it], etime = stime[it+1])


os.system(call_matlab + "\"ElSpec_IC_iter(" + str(i) + ", \\\"" + log_directory + "\\\", \\\"" + path_eiscat_data + "\\\")\" -nodisplay")

while True:
    ic4elspec.ic(log_directory, chemistry_config, 'ElSpec-iqt_IC_', i, mixf = mixf)
    i = i+1
    print('Mean Iteration Duration:', (time.time() - start_t)/i, 's')
    os.system(call_matlab + "\"ElSpec_IC_iter(" + str(i) + ", \\\"" + log_directory + "\\\", \\\"" + path_eiscat_data + "\\\")\" -nodisplay")


"""
to adapt this script for processing batches of 15minutes:
 1. divide data set in batches
 2. adapt this script for batches
 2. adapt ic4elspec such that
    a. time starts at the proper time (no 30 min beforehand, or only optionally)
    b. 
3. adapt matlab function such that time intervals can be set, and initial parameters are loaded (ie ne, composition)
"""


"""


import matlab.engine
m = matlab.engine.start_matlab()
m.addpath("../ELSPEC")

import numpy as np
egrid = np.logspace(1,5,200)
fitdir = '../Data/Eiscat/fit'
ppdir = '../Data/Eiscat/pp'
experiment = 'arc1'
hmax = 150
hmin = 95
btime = np.array([2006, 12, 12, 19, 30, 0], dtype=float)
etime = np.array([2006, 12, 12, 19, 30, 10], dtype=float)
ionomodel = 'Sergienko'
recombmodel = 'SheehanGr'
#recombmodel = ['SheehanGrFlipchem']
integtype = 'integrate'
tres = 'best'
maxorder = 5
ninteg = 20
ErrType = 'l'
Outname = m.fullfile(log_directory, ["ElSpec-iqt_IC_" + str(i)])[0]
print(Outname)


m.ElSpec_iqt_ic('fitdir',fitdir,
                'ppdir', ppdir,
                'experiment', experiment,
                'hmax', hmax, 'hmin', hmin,
                'btime', btime, 'etime', etime,
                'ionomodel', ionomodel,
                'integtype', integtype,
                'egrid', egrid,
                'tres', tres,
                'ErrType', ErrType,
                'MaxOrder', maxorder,
                'ninteg', ninteg,
                'Outfilename', Outname,
                'iri_ic', 0,
                'alpha_eff', 0,
                'iteration', i,
                'ne_init', 0,
                'recombmodel', recombmodel,
                nargout = 0
                )
"""
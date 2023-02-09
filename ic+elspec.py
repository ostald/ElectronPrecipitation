import os
import import_ipynb
import ic4elspec
import time
#%matplotlib widget

i = 0
#os.system("/Applications/MATLAB_R2022b.app/bin/./matlab -sd"+
#          " \"/Users/ost051/Documents/PhD/ELSPEC-2022" +
#          "\" -batch \"ElSpec_IC\" -nodisplay")

ic4elspec.ic('/Users/ost051/Documents/PhD/ELSPEC-2022/Reactions Set Rees/',
             'ElSpec-iqt_IC_', i)
i = 1

while True:
    osstr = "/Applications/MATLAB_R2022b.app/bin/./matlab -sd "+\
            " \"/Users/ost051/Documents/PhD/ELSPEC-2022"+\
            "\" -batch \"ElSpec_IC_iter("+str(i)+")\" -nodisplay"
    os.system(osstr)
    ic4elspec.ic('/Users/ost051/Documents/PhD/ELSPEC-2022/Reactions Set Rees/',
                 'ElSpec-iqt_IC_', i)
    i = i+1 




import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat
import ionChem

f1 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.02.16_17_58_37 mixf=0/IC_res_8.pickle'
f2 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.02.16_21_45_58 mixf=1/IC_res_66.pickle'

with open(f1, 'rb') as pf1:
    [ts1, z1, n_ic1, eff_rr1] = pickle.load(pf1)

with open(f2, 'rb') as pf2:
    [ts2, z2, n_ic2, eff_rr2] = pickle.load(pf2)

d_n_ic = n_ic1 - n_ic2
d_eff_rr = eff_rr1 - eff_rr2

fig, axs = plt.subplots()
pc = axs.pcolormesh(ts1, z1, d_eff_rr, norm=mpl.colors.CenteredNorm(),
                          label='alpha', cmap='RdBu')
fig.suptitle('Difference in Rec. Rate')
fig.supylabel('Altitude [km]')
fig.supxlabel('Time [s]')
plt.colorbar(pc, ax = axs)
plt.show()



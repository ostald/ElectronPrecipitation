import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat
import ionChem

f1 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.02.16_17_58_37 mixf=0/IC_res_8.pickle'
f2 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.02.20_17_18_54 mixf=1/IC_res_49.pickle'

f1 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.02.20_17_16_29 mixf=1/IC_res_0.pickle'
f2 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.02.20_17_16_29 mixf=1/IC_res_48.pickle'

spstr = ['e','O','O2','O2p','Np','N2','N2p','NO','NOp','H','Hp','O_1D','O_1S','N_2D','N_4S','O2p_a4P','Op_2D', \
         'Op_4S','Op_2P']

with open(f1, 'rb') as pf1:
    [ts1, z1, n_ic1, eff_rr1] = pickle.load(pf1)

with open(f2, 'rb') as pf2:
    [ts2, z2, n_ic2, eff_rr2] = pickle.load(pf2)

d_n_ic = n_ic1 - n_ic2
d_eff_rr = eff_rr1 - eff_rr2

fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
pc0 = axs.flat[0].pcolormesh(ts1, z1, eff_rr1)#, norm=mpl.colors.CenteredNorm(),
#                          label='alpha', cmap='RdBu')
pc1 = axs.flat[1].pcolormesh(ts2, z2, eff_rr2)#, norm=mpl.colors.CenteredNorm(),
#                          label='alpha', cmap='RdBu')
pc2 = axs.flat[2].pcolormesh(ts1, z1, d_eff_rr, norm=mpl.colors.CenteredNorm(),
                          label='alpha', cmap='RdBu')
fig.suptitle('Difference in Rec. Rate')
fig.supylabel('Altitude [km]')
fig.supxlabel('Time [s]')
plt.colorbar(pc0, ax = axs.flat[0])
plt.colorbar(pc1, ax = axs.flat[1])
plt.colorbar(pc2, ax = axs.flat[2])

for i in range(15):
    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
    pc0 = axs.flat[0].pcolormesh(ts1, z1, n_ic1[:, i, :], norm=mpl.colors.LogNorm())
    pc1 = axs.flat[1].pcolormesh(ts2, z2, n_ic2[:, i, :], norm=mpl.colors.LogNorm())
    pc2 = axs.flat[2].pcolormesh(ts1, z1, d_n_ic[:, i, :]/n_ic1[:, i, :], norm=mpl.colors.CenteredNorm(),
                                 label='alpha', cmap='RdBu')
    axs.flat[0].set_title('Rel. Difference in ' + spstr[i])
    fig.supylabel('Altitude [km]')
    fig.supxlabel('Time [s]')
    plt.colorbar(pc0, ax = axs.flat[0])
    plt.colorbar(pc1, ax = axs.flat[1])
    plt.colorbar(pc2, ax = axs.flat[2])
plt.show()



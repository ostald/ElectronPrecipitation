import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat

plt.rcParams.update({'font.size': 10})

f = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/ElSpec-iqt_IC_7.mat'
con = loadmat.loadmat(f)["ElSpecOut"]

ne_meas = con['pp']
ne_std = con['ppstd']
ne_mod = con['ne']
ts = con['ts']
h = con['h']

ts = ts - ts[0]

res = (ne_meas - ne_mod) / ne_std

fig, axs = plt.subplots(figsize=(6.4, 3))
pc0 = axs.pcolormesh(ts, h, res, norm = mpl.colors.CenteredNorm(halfrange = 5), cmap = plt.cm.seismic)
plt.colorbar(pc0, ax=axs, label='Norm. Residuals [1]', anchor=(-0.25,0.0))
fig.supylabel('Altitude [km]')
fig.supxlabel('Time [s]')
plt.title('Normalized Residuals in $\mathrm{n_e}$')
plt.tight_layout()
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/residuals.png')
plt.show()

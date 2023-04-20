import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat
import ionChem

direc = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.04.18_14_33_47 mixf=0/'
files = glob.glob(direc + '*.pickle')

elspec_0 = loadmat.loadmat(direc + 'ElSpec-iqt_IC_0.mat')["ElSpecOut"]

nax = 13

all_data = []

for i, f in enumerate(files[:nax]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    print(f)
    with open(f, 'rb') as pf:
        data = pickle.load(pf)
        all_data = [*all_data, data[:4]]

fig,ax = plt.subplots()
for i, f in enumerate(files[:nax]):
    if True:
        data = all_data[i]
        eff_rr = data[3][:, 1:]
        #eff_rr = data[2][:, 0, 1:]
    if i>0:
        if True:
            data = all_data[i-1]
            eff_rr_o = data[3][:, 1:]
            #eff_rr_o = data[2][:, 0, 1:]
    else:
        eff_rr_o = elspec_0["alpha"]
        #eff_rr_o = elspec_0["ne"]
    d_effrr_r = (eff_rr_o - eff_rr) / eff_rr
    if i == 0:
        ax.plot(i, np.sum(np.abs(d_effrr_r.flat))/len(d_effrr_r.flat), 'x', color = 'black', label = 'Mean')
        ax.plot(i, np.max(np.abs(d_effrr_r.flat)), 'x', color = 'blue', label = 'Max')
    else:
        ax.plot(i, np.sum(np.abs(d_effrr_r.flat)) / len(d_effrr_r.flat), 'x', color='black')
        ax.plot(i, np.max(np.abs(d_effrr_r.flat)), 'x', color='blue')
plt.title(r' Relative Variation in $\alpha_{eff}$')
#plt.title(r' Relative Variation in $n_e$')
ax.set_ylabel(r"$\frac{\alpha_{eff, i-1}}{\alpha_{eff, i}} - 1$")
#ax.set_ylabel(r"$\frac{n_{e, i-1}}{n_{e, i}} - 1$")
ax.set_xlabel("Iteration i")
ax.set_yscale('log')
ax.legend()
plt.savefig('/Users/ost051/Documents/PhD/Electron Precipitation/writing/plots/alpha_rel_dev_mixf0.png')
#plt.savefig('/Users/ost051/Documents/PhD/Electron Precipitation/writing/plots/ne_rel_dev.png')
plt.show()

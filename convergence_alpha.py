import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat
import mat73

#direc = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/'
direc = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.10.03_12_05_57testing_convergence_lsqnonlin_mixf=0/'
#files = glob.glob(direc + '*.pickle')
files = glob.glob(direc + 'IC*')

try:
    elspec_0 = loadmat.loadmat(direc + 'ElSpec-iqt_IC_0.mat')["ElSpecOut"]
except:
    elspec_0 = mat73.loadmat(direc + 'ElSpec-iqt_IC_0.mat')["ElSpecOut"]

nax = 10

all_data = []

plt.rcParams.update({'font.size': 12})

# for i, f in enumerate(files[:nax]):
#     f = direc + 'IC_res_'+str(i)+'.pickle'
#     print(f)
#     with open(f, 'rb') as pf:
#         data = pickle.load(pf)
#         all_data = [*all_data, data[:4]]

for i, f in enumerate(files[:nax]):
    f = direc + 'IC_'+str(i)+'.mat'
    print(f)
    try:
        data = loadmat.loadmat(f)
    except:
        data = mat73.loadmat(f)
    all_data = [*all_data, data]

fig,ax = plt.subplots(figsize=(6.4, 5))
for i, f in enumerate(files[:nax]):
    data = all_data[i]
    #eff_rr = data[3][:, 1:]
    eff_rr = data['eff_rr']
    if i>0:
        data = all_data[i-1]
        #eff_rr_o = data[3][:, 1:]
        eff_rr_o = data['eff_rr']
    else:
        eff_rr_o = elspec_0["alpha"]
    d_effrr_r = (eff_rr_o - eff_rr) / eff_rr
    if i == 0:
        ax.plot(i+1, np.sum(np.abs(d_effrr_r.flat))/len(d_effrr_r.flat), 'x', color = 'black', label = 'Mean')
        ax.plot(i+1, np.max(np.abs(d_effrr_r.flat)), 'x', color = 'blue', label = 'Max')
    else:
        ax.plot(i+1, np.sum(np.abs(d_effrr_r.flat)) / len(d_effrr_r.flat), 'x', color='black')
        ax.plot(i+1, np.max(np.abs(d_effrr_r.flat)), 'x', color='blue')
    ax.plot(0, 0, 'x', alpha=0)
plt.title(r'Relative Variation in $\alpha_{eff}$ between Iterations')
ax.set_ylabel(r"$(\alpha_{eff, i-1}/\alpha_{eff, i}) - 1$")
ax.set_xlabel("Iteration i")
ax.set_yscale('log')
ax.legend()
#plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/alpha_rel_dev_mixf0.png')
plt.show()
exit()

fig,ax = plt.subplots(figsize=(6.4, 5))
for i, f in enumerate(files[:nax]):
    data = all_data[i]
    ne = data[2][:, 0, 1:]
    if i>0:
        data = all_data[i-1]
        ne_o = data[2][:, 0, 1:]
    else:
        ne_o = elspec_0["ne"]
    d_effrr_r = (ne_o - ne) / ne
    if i == 0:
        ax.plot(i+1, np.sum(np.abs(d_effrr_r.flat))/len(d_effrr_r.flat), 'x', color = 'black', label = 'Mean')
        ax.plot(i+1, np.max(np.abs(d_effrr_r.flat)), 'x', color = 'blue', label = 'Max')
    else:
        ax.plot(i+1, np.sum(np.abs(d_effrr_r.flat)) / len(d_effrr_r.flat), 'x', color='black')
        ax.plot(i+1, np.max(np.abs(d_effrr_r.flat)), 'x', color='blue')
    ax.plot(0, 0, 'x', alpha = 0)
plt.title(r' Relative Variation in Electron Density between Iterations')
ax.set_ylabel(r"$(n_{e, i-1}/n_{e, i}) - 1$")
ax.set_xlabel("Iteration i")
ax.set_yscale('log')
ax.legend()
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/ne_rel_dev.png')
plt.show()

import convergence_bjorn
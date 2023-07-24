import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import loadmat

mpl.rc('image', cmap = 'plasma')

mat0 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/ElSpec-iqt_IC_0.mat'
f0 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/IC_res_0.pickle'


mat2 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/ElSpec-iqt_IC_14.mat'
f2 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/IC_res_14.pickle'


con = loadmat.loadmat(mat0)["ElSpecOut"]
ts0 = con['ts']
ts0 = ts0 - ts0[0]
z0 = con["h"]


con2 = loadmat.loadmat(mat2)["ElSpecOut"]

def plot_compare(x, y, data1, data2, title, label, vminlimit = None):
    vmin = min(data1.min(), data2.min())
    vmax = max(data1.max(), data2.max())

    if vminlimit is not None:
        if vminlimit > vmin:
            vmin = vminlimit
            data1[data1 < vminlimit] = 0
            data2[data2 < vminlimit] = np.nan
    print(vmin, vmax)

    diff = data2 - data1
    rdiff = diff / data2

    norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize =(6.4, 5))
    pc0 = axs.flat[0].pcolormesh(x, y, data1, norm=norm)
    pc1 = axs.flat[1].pcolormesh(x, y, data2, norm=norm)
    pc2 = axs.flat[2].pcolormesh(x, y, rdiff, norm=mpl.colors.CenteredNorm(), \
                             label='alpha', cmap=mpl.cm.RdBu.reversed())
    #pc3 = axs.flat[3].pcolormesh(x, y, n_ic2[:, 0, 1:], norm = mpl.colors.LogNorm())
    #pc4 = axs.flat[4].pcolormesh(x, y, q_final, norm = mpl.colors.LogNorm())

    #fig.suptitle(title)
    fig.supylabel('Altitude [km]')
    fig.supxlabel('Time [s]')
    axs.flat[0].text(5, 142, 'a')
    axs.flat[1].text(5, 142, 'b')
    axs.flat[2].text(5, 142, 'c')
    #axs.flat[0].set_title('Initial Values')
    #axs.flat[1].set_title('Final Values')
    #axs.flat[2].set_title('Relative Difference')
    plt.tight_layout()
    try:
        plt.colorbar(pc0, ax = axs.flat[0:2], label = label, anchor=(-0.25,0.0))
        plt.colorbar(pc2, ax = axs.flat[2], label = 'Relative Variation [1]', anchor=(-0.25,0.0))
        #plt.colorbar(pc3, ax = axs.flat[3], label = 'ne', anchor=(-0.25,0.0))
        #plt.colorbar(pc4, ax = axs.flat[4], label = r'$q_{e}$', anchor=(-0.25,0.0))
    except ValueError:
        print(vmin, vmax)

plot_compare(ts0, z0, con["ne"], con2["ne"], 'e- Density', r'$\mathrm{n_{e} \, [m^{-3}]}$')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/ne_elspec_0_end.png')

plot_compare(ts0, z0, con["alpha"], con2["alpha"], 'Effective Recombination Rate', r'$\mathrm{\alpha_{eff} \, [m^{3}s^{-1}]}$')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/comp_alpha_2.png')


def plot_compare_alpha(x, y, data1, data2, title, label, vminlimit = None):
    vmin = min(data1.min(), data2.min())
    vmax = max(data1.max(), data2.max())

    vmine = min(con2['ne'].min(), con2['pp'].min())
    vmaxe = max(con2['ne'].max(), con2['pp'].max())
    if vmine < 1e7: vmine = 1e9

    if vminlimit is not None:
        if vminlimit > vmin:
            vmin = vminlimit
            data1[data1 < vminlimit] = 0
            data2[data2 < vminlimit] = np.nan
    print(vmin, vmax)

    diff = data2 - data1
    rdiff = diff / data2

    norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    norme = mpl.colors.LogNorm(vmin=vmine, vmax=vmaxe)

    fig, axs = plt.subplots(nrows=6, ncols=1, sharex=True, figsize =(6.4, 10))
    pc0 = axs.flat[0].pcolormesh(x, y, data1, norm=norm)
    pc1 = axs.flat[1].pcolormesh(x, y, data2, norm=norm)
    pc2 = axs.flat[2].pcolormesh(x, y, rdiff, norm=mpl.colors.CenteredNorm(), \
                             label='alpha', cmap=mpl.cm.RdBu.reversed())
    pc3 = axs.flat[3].pcolormesh(x, y, con2['ne'], norm=norme)#mpl.colors.LogNorm())
    pc4 = axs.flat[4].pcolormesh(x, y, con2['pp'], norm=norme)#mpl.colors.LogNorm())
    pc5 = axs.flat[5].pcolormesh(x, y, con2['par'][:, 2, :])  # , norm = mpl.colors.LogNorm())

    #fig.suptitle(title)
    fig.supylabel('Altitude [km]')
    fig.supxlabel('Time [s]')
    axs.flat[0].text(5, 142, 'a')
    axs.flat[1].text(5, 142, 'b')
    axs.flat[2].text(5, 142, 'c')
    axs.flat[3].text(5, 142, 'd')
    axs.flat[4].text(5, 142, 'e')
    axs.flat[5].text(5, 142, 'f')
    #axs.flat[0].set_title('Initial Values')
    #axs.flat[1].set_title('Final Values')
    #axs.flat[2].set_title('Relative Difference')
    plt.tight_layout()
    try:
        plt.colorbar(pc0, ax = axs.flat[0:2], label = label, anchor=(-0.25,0.0))
        plt.colorbar(pc2, ax = axs.flat[2], label = 'Relative Variation [1]', anchor=(-0.25,0.0))
        plt.colorbar(pc3, ax=axs.flat[3], label='$\mathrm{n_{e, mod} [m^{-3}]}$', anchor=(-0.25,0.0))
        plt.colorbar(pc4, ax=axs.flat[4], label='$\mathrm{n_{e, meas} [m^{-3}]}$', anchor=(-0.25,0.0))
        plt.colorbar(pc5, ax=axs.flat[5], label=r'$\mathrm{T_{e} [K]}$', anchor=(-0.25,0.0))
    except ValueError:
        print(vmin, vmax)
plot_compare_alpha(ts0, z0, con["alpha"], con2["alpha"], 'Effective Recombination Rate', r'$\mathrm{\alpha_{eff} \, [m^{3}s^{-1}]}$')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/comp_alpha.png')

n_model = con["iri"]
[Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
[nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * con["ne"]
with open(f2, 'rb') as pf2:
    data2 = pickle.load(pf2)
n_ic2 = data2[2]
plot_compare(ts0, z0, nO2p/nNOp, n_ic2[:, 3, 1:]/n_ic2[:, 8, 1:], 'O2+/NO+ Density Ratio', r'$\mathrm{n(O_2^+)/n(NO^+)} \, [1]$')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/ratio_o2+_no+.png')


def plot_compare_ionratio(x, y, data1, data2, title, label, vminlimit = None):
    vmin = min(data1.min(), data2.min())
    vmax = max(data1.max(), data2.max())

    if vminlimit is not None:
        if vminlimit > vmin:
            vmin = vminlimit
            data1[data1 < vminlimit] = 0
            data2[data2 < vminlimit] = np.nan
    print(vmin, vmax)

    diff = data2 - data1
    rdiff = diff / data2

    norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize =(6.4, 4))
    pc0 = axs.flat[0].pcolormesh(x, y, data1, norm=norm)
    pc1 = axs.flat[1].pcolormesh(x, y, data2, norm=norm)
    #pc2 = axs.flat[2].pcolormesh(x, y, rdiff, norm=mpl.colors.CenteredNorm(), \
    #                         label='alpha', cmap=mpl.cm.RdBu.reversed())
    #pc3 = axs.flat[3].pcolormesh(x, y, n_ic2[:, 0, 1:], norm = mpl.colors.LogNorm())
    #pc4 = axs.flat[4].pcolormesh(x, y, q_final, norm = mpl.colors.LogNorm())

    #fig.suptitle(title)
    fig.supylabel('Altitude [km]')
    fig.supxlabel('Time [s]')
    axs.flat[0].text(5, 142, 'a')
    axs.flat[1].text(5, 142, 'b')
    #axs.flat[2].text(5, 142, 'c')
    #axs.flat[0].set_title('Initial Values')
    #axs.flat[1].set_title('Final Values')
    #axs.flat[2].set_title('Relative Difference')
    plt.tight_layout()
    try:
        plt.colorbar(pc0, ax = axs.flat[0:2], label = label, anchor=(-0.25,0.0))
        #plt.colorbar(pc2, ax = axs.flat[2], label = 'Relative Variation [1]', anchor=(-0.25,0.0))
        #plt.colorbar(pc3, ax = axs.flat[3], label = 'ne', anchor=(-0.25,0.0))
        #plt.colorbar(pc4, ax = axs.flat[4], label = r'$q_{e}$', anchor=(-0.25,0.0))
    except ValueError:
        print(vmin, vmax)

n_model = con["iri"]
[Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
[nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * con["ne"]
with open(f2, 'rb') as pf2:
    data2 = pickle.load(pf2)
n_ic2 = data2[2]
plot_compare_ionratio(ts0, z0, nO2p/nNOp, n_ic2[:, 3, 1:]/n_ic2[:, 8, 1:], 'O2+/NO+ Density Ratio', r'$\mathrm{n(O_2^+)/n(NO^+)} \, [1]$')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/ratio_o2+_no+_v2.png')



mat3 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/ElSpec-iqt_IC_1.mat'
con3 = loadmat.loadmat(mat3)["ElSpecOut"]
f3 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/IC_res_0.pickle'
n_model = con["iri"]
[Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
[nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * con["ne"]
with open(f3, 'rb') as pf2:
    data2 = pickle.load(pf2)
n_ic2 = data2[2]
ts2 = data2[0]
plot_compare(ts2[1:], z0, nO2p/nNOp, n_ic2[:, 3, 1:]/n_ic2[:, 8, 1:], 'O2+/NO+ Density Ratio (first iteration)', r'$\mathrm{n(O_2^+)/n(NO^+)} \, [1]$', vminlimit=0.1)
#plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/ratio_o2+_no+.png')
#plt.show()


with open(f2, 'rb') as pf2:
    data2 = pickle.load(pf2)
    [ts2, z2, n_ic2, eff_rr2] = data2[:4]


# norm = mpl.colors.LogNorm()
fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(6.4, 5))
pc0 = axs.flat[0].pcolormesh(ts2[1:], z2, n_ic2[:, 3, 1:], norm=mpl.colors.LogNorm())
pc1 = axs.flat[1].pcolormesh(ts2[1:], z2, n_ic2[:, 8, 1:], norm=mpl.colors.LogNorm())
pc2 = axs.flat[2].pcolormesh(ts2[1:], z2, n_ic2[:,17, 1:], norm=mpl.colors.LogNorm())
plt.colorbar(pc0, ax=axs.flat[0], label=r'$\mathrm{n_{O_2^+} [m^{-3}]}$', anchor=(-0.25,0.0))
plt.colorbar(pc1, ax=axs.flat[1], label=r'$\mathrm{n_{NO^+} [m^{-3}]}$', anchor=(-0.25,0.0))
plt.colorbar(pc2, ax=axs.flat[2], label=r'$\mathrm{n_{O^+} [m^{-3}]}$', anchor=(-0.25,0.0))
fig.supylabel('Altitude [km]')
fig.supxlabel('Time [s]')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/n_ions.png')
#plt.show()



fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(6.4, 5))
pc0 = axs.flat[0].pcolormesh(ts2[1:], z2, n_ic2[:, 3, 1:]/n_ic2[:, 0, 1:], vmin=0, vmax=1)#, norm=mpl.colors.LogNorm())
pc1 = axs.flat[1].pcolormesh(ts2[1:], z2, n_ic2[:, 8, 1:]/n_ic2[:, 0, 1:], vmin=0, vmax=1)#, norm=mpl.colors.LogNorm())
pc2 = axs.flat[2].pcolormesh(ts2[1:], z2, n_ic2[:,17, 1:]/n_ic2[:, 0, 1:], vmin=0, vmax=1)#, norm=mpl.colors.LogNorm())
plt.colorbar(pc0, ax=axs.flat[0], label=r'$\mathrm{n_{O_2^+}/n_e \, [1]}$', anchor=(-0.25,0.0))
plt.colorbar(pc1, ax=axs.flat[1], label=r'$\mathrm{n_{NO^+}/n_e \, [1]}$', anchor=(-0.25,0.0))
plt.colorbar(pc2, ax=axs.flat[2], label=r'$\mathrm{n_{O^+}/n_e \, [1]}$', anchor=(-0.25,0.0))
fig.supylabel('Altitude [km]')
fig.supxlabel('Time [s]')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/n_ions_rel_i13.png')


fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(6.4, 5))
pc0 = axs.flat[0].pcolormesh(ts0, z0, con["iri"][:, 8, :]/np.sum(con["iri"][:, 7:, :], axis = 1), vmin=0, vmax=1)#, norm=mpl.colors.LogNorm())
pc1 = axs.flat[1].pcolormesh(ts0, z0, con["iri"][:, 7, :]/np.sum(con["iri"][:, 7:, :], axis = 1), vmin=0, vmax=1)#, norm=mpl.colors.LogNorm())
pc2 = axs.flat[2].pcolormesh(ts0, z0, con["iri"][:, 9, :]/np.sum(con["iri"][:, 7:, :], axis = 1), vmin=0, vmax=1)#, norm=mpl.colors.LogNorm())
plt.colorbar(pc0, ax=axs.flat[0], label=r'$\mathrm{n_{O_2^+}/n_e \, [1]}$', anchor=(-0.25,0.0))
plt.colorbar(pc1, ax=axs.flat[1], label=r'$\mathrm{n_{NO^+}/n_e \, [1]}$', anchor=(-0.25,0.0))
plt.colorbar(pc2, ax=axs.flat[2], label=r'$\mathrm{n_{O^+}/n_e \, [1]}$', anchor=(-0.25,0.0))
fig.supylabel('Altitude [km]')
fig.supxlabel('Time [s]')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/n_ions_rel_i0.png')



fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(6.4, 5))

egrid = con["egrid"]

Ie0 = con["Ie"]
pc0 = axs.flat[0].pcolormesh(ts0, egrid[:-1] / 1e3, Ie0*np.tile(egrid[:-1], (674, 1)).T, norm= mpl.colors.LogNorm(vmin = 1e8, vmax = 1e13))
Ie2 = con2["Ie"]
pc1 = axs.flat[1].pcolormesh(ts0, egrid[:-1] / 1e3, Ie2*np.tile(egrid[:-1], (674, 1)).T, norm= mpl.colors.LogNorm(vmin = 1e8, vmax = 1e13))
pc2 = axs.flat[2].pcolormesh(ts0, egrid[:-1] / 1e3, (Ie2 - Ie0)*np.tile(egrid[:-1], (674, 1)).T, norm = mpl.colors.CenteredNorm(halfrange=5e12), cmap = plt.cm.RdBu.reversed())
for ax in axs:
    ax.set_yscale('log')
    ax.set_ylim(1, 100)

plt.colorbar(pc0, ax=axs.flat[0:2], label = r'$\mathrm{\phi * E \, [eV eV^{-1} m^{-2} s^{-1}]}$', anchor=(-0.25,0.0))
#plt.colorbar(pc1, ax=axs.flat[1], label = r'$\mathrm{\phi * E \, [eV eV^{-1} m^{-2} s^{-1}]}$', anchor=(-0.25,0.0))
plt.colorbar(pc2, ax=axs.flat[2], label = r'$\mathrm{\phi * E \, [eV eV^{-1} m^{-2} s^{-1}]}$', anchor=(-0.25,0.0))
fig.supylabel('Energy [eV]')
fig.supxlabel('Time [s]')

plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/compare_energies.png')
plt.show()



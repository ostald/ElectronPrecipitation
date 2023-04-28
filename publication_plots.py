import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import loadmat

mat0 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.04.27_17_54_27_mixf=0/ElSpec-iqt_IC_0.mat'

mat2 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.04.27_17_54_27_mixf=0/ElSpec-iqt_IC_14.mat'
f2 = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.04.27_17_54_27_mixf=0/IC_res_14.pickle'


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
                             label='alpha', cmap='RdBu')
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
        plt.colorbar(pc0, ax = axs.flat[0:2], label = label)
        plt.colorbar(pc2, ax = axs.flat[2], label = [1])
        #plt.colorbar(pc3, ax = axs.flat[3], label = 'ne')
        #plt.colorbar(pc4, ax = axs.flat[4], label = r'$q_{e}$')
    except ValueError:
        print(vmin, vmax)

plot_compare(ts0, z0, con["ne"], con2["ne"], 'e- Density', r'$\mathrm{n_{e} \, [m^{3}s^{-1}]}$')
plt.savefig('/Users/ost051/Documents/PhD/Electron Precipitation/writing/plots/ne_elspec_0_end.png')


n_model = con["iri"]
[Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
[nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * con2["ne"]
with open(f2, 'rb') as pf2:
    data2 = pickle.load(pf2)
n_ic2 = data2[2]
plot_compare(ts0, z0, nO2p/nNOp, n_ic2[:, 3, 1:]/n_ic2[:, 8, 1:], 'O2+/NO+ Density Ratio', r'$\mathrm{n(O_2^+)/n(NO^+)} \, [1]$')
plt.show()
plt.savefig('/Users/ost051/Documents/PhD/Electron Precipitation/writing/plots/ratio_o2+_no+.png')



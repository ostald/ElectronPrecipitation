import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat


eval_list = [
'2023.04.28_09_33_05_mixf=1',
'2023.05.03_14_07_04_mixf=0',
'2023.05.08_18_31_02_mixf=0',
'2023.05.10_17_18_42_mixf=0',
'2023.05.12_09_56_23_mixf=0'
]

for dir in eval_list:

    f1 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/IC_res_7.pickle'
    mat0 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/ElSpec-iqt_IC_7.mat'

    f2 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/'+ dir +'/IC_res_7.pickle'
    mat2 = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/'+ dir +'/ElSpec-iqt_IC_7.mat'

    spstr = ['e','O','O2','O2p','Np','N2','N2p','NO','NOp','H','Hp','O_1D','O_1S','N_2D','N_4S','O2p_a4P','Op_2D', \
             'Op_4S','Op_2P']


    con = loadmat.loadmat(mat0)["ElSpecOut"]
    n_model = con["iri"]
    [Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
    eff_rr0 = con['alpha']
    ts0 = con['ts']
    ts0 = ts0 - ts0[0]
    z0 = con["h"]
    ne = con["ne"]
    q0 = con["q"]
    Ie0 = con["Ie"]
    E = con["E"][:-1]
    [nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * ne
    [ne_, Ti, Te, _] = con["par"].swapaxes(0, 1)

    con2 = loadmat.loadmat(mat2)["ElSpecOut"]
    q_final = con2["q"]
    Ief = con2["Ie"]



    with open(f1, 'rb') as pf1:
        data1 = pickle.load(pf1)
        [ts1, z1, n_ic1, eff_rr1] = data1[:4]

    with open(f2, 'rb') as pf2:
        data2 = pickle.load(pf2)
        [ts2, z2, n_ic2, eff_rr2] = data2[:4]

    dnNOp = nNOp - n_ic2[:, 8, 1:]
    dnO2p = nO2p - n_ic2[:, 3, 1:]
    dnOp  = nOp  - n_ic2[:, 16, 1:]

    #norm = mpl.colors.LogNorm()
    # fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize =(6.4, 5))
    # pc0 = axs.flat[0].pcolormesh(ts2[1:], z2, n_ic2[:, 3, 1:], norm=mpl.colors.LogNorm())
    # pc1 = axs.flat[1].pcolormesh(ts2[1:], z2, n_ic2[:, 8, 1:], norm=mpl.colors.LogNorm())
    # pc2 = axs.flat[2].pcolormesh(ts2[1:], z2, n_ic2[:,17, 1:], norm=mpl.colors.LogNorm())
    # plt.colorbar(pc0, ax = axs.flat[0], label=r'$\mathrm{n_{O_2^+ [m^{-3}]}}$')
    # plt.colorbar(pc1, ax = axs.flat[1], label=r'$\mathrm{n_{NO^+ [m^{-3}]}}$')
    # plt.colorbar(pc2, ax = axs.flat[2], label=r'$\mathrm{n_{O^+ [m^{-3}]}}$')
    # #fig.suptitle('asdf')
    # fig.supylabel('Altitude [km]')
    # fig.supxlabel('Time [s]')
    #plt.show()
    #plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/n_ions.png')



    Ie0 = np.array([i*e for i, e in zip(Ie0, E)])
    Ief = np.array([i*e for i, e in zip(Ief, E)])

    vmin = min(Ie0.min(), Ief.min())
    vmax = max(Ie0.max(), Ief.max())
    vminlimit = 1
    if vminlimit is not None:
        if vminlimit > vmin:
            vmin = vminlimit
            Ie0[Ie0 < vminlimit] = vminlimit
            Ief[Ief < vminlimit] = vminlimit
    print(vmin, vmax)

    diff = Ief - Ie0
    rdiff = diff / Ief

    norm = mpl.colors.LogNorm(vmin=1e7, vmax=vmax)

    # fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True, figsize =(6.4, 8))
    #
    # plt.gcf().axes[0].set_yscale('log')
    # plt.gcf().axes[1].set_yscale('log')
    # #plt.gcf().axes[0].set_ylim([1, 1e6])
    # #plt.gcf().axes[1].set_ylim([1, 1e6])
    #
    # pc0 = axs.flat[0].pcolormesh(ts0, E, Ie0, norm=norm)
    # pc1 = axs.flat[1].pcolormesh(ts0, E, Ief, norm=norm)
    # #pc2 = axs.flat[2].pcolormesh(ts0, E, rdiff, norm=mpl.colors.CenteredNorm(), \
    # #                             label='alpha', cmap='RdBu')
    # pc3 = axs.flat[3].pcolormesh(ts0, z0, n_ic2[:, 0, 1:], norm = mpl.colors.LogNorm())
    # pc4 = axs.flat[4].pcolormesh(ts0, z0, q_final, norm = mpl.colors.LogNorm())
    #
    # #fig.suptitle(title)
    # axs.flat[0].set_ylabel('Energy [eV]')
    # axs.flat[1].set_ylabel('Energy [eV]')
    # axs.flat[2].set_ylabel('Altitude [km]')
    # axs.flat[3].set_ylabel('Altitude [km]')
    # axs.flat[4].set_ylabel('Altitude [km]')
    # fig.supxlabel('Time [s]')
    # axs.flat[0].set_title('Initial Values')
    # axs.flat[1].set_title('Final Values')
    # axs.flat[2].set_title('Relative Difference')
    # plt.tight_layout()
    # try:
    #     plt.colorbar(pc0, ax = axs.flat[0:2], label = 'a')
    #     #plt.colorbar(pc2, ax = axs.flat[2], label = [1])
    #     plt.colorbar(pc3, ax = axs.flat[3], label = 'ne')
    #     plt.colorbar(pc4, ax = axs.flat[4], label = r'$q_{e}$')
    # except ValueError:
    #     print(vmin, vmax)
    #
    # plt.show()

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

        fig, axs = plt.subplots(nrows=6, ncols=1, sharex=True, figsize =(6.4, 8.5))
        fig.canvas.manager.set_window_title(dir)
        pc0 = axs.flat[0].pcolormesh(x, y, data1, norm=norm)
        pc1 = axs.flat[1].pcolormesh(x, y, data2, norm=norm)
        pc2 = axs.flat[2].pcolormesh(x, y, rdiff, norm=mpl.colors.CenteredNorm(), \
                                     label='alpha', cmap='RdBu')
        pc3 = axs.flat[3].pcolormesh(x, y, n_ic2[:, 0, 1:], norm = mpl.colors.LogNorm())
        pc4 = axs.flat[4].pcolormesh(x, y, q_final, norm = mpl.colors.LogNorm())
        pc5 = axs.flat[5].pcolormesh(x, y, Te)#, norm = mpl.colors.LogNorm())

        fig.suptitle(title)
        fig.supylabel('Altitude [km]')
        fig.supxlabel('Time [s]')
        axs.flat[0].text(5, 144, 'a')
        axs.flat[1].text(5, 144, 'b')
        axs.flat[3].text(5, 144, 'c')
        #axs.flat[0].set_title('Initial Values')
        #axs.flat[1].set_title('Final Values')
        #axs.flat[2].set_title('Relative Difference')
        plt.tight_layout()
        try:
            plt.colorbar(pc0, ax = axs.flat[0:2], label = label)
            plt.colorbar(pc2, ax = axs.flat[2], label = [1])
            plt.colorbar(pc3, ax = axs.flat[3], label = 'ne')
            plt.colorbar(pc4, ax = axs.flat[4], label = r'$q_{e}$')
            plt.colorbar(pc5, ax = axs.flat[5], label = r'$T_{e}$')
        except ValueError:
            print(vmin, vmax)


    #plot_compare(ts0, z0, q0, q_final, 'Production Rate', r'$q_{ne} \, [m^{3}s^{-1}$]')

    plot_compare(ts0, z0, eff_rr0, eff_rr2[:, 1:], 'Recombination Rate', r'$\mathrm{\alpha_{eff} \, [m^{-3}s^{-1}]}$')

    plot_compare(ts0, z0, con["ne"], con2["ne"], 'e- Density', r'$\mathrm{n(e^-) \, [m^{3}s^{-1}]}$')

    plot_compare(ts0, z0, nNOp, n_ic2[:, 8, 1:], 'NO+ Density', r'$\mathrm{n(NO^+) \, [m^{3}s^{-1}]}$')

    plot_compare(ts0, z0, nO2p, n_ic2[:, 3, 1:], 'O2+ Density', r'$\mathrm{n(O_2^+) \, [m^{3}s^{-1}]}$')

    plot_compare(ts0, z0, nOp, n_ic2[:, 17, 1:], 'O+ Density', r'$\mathrm{n(O^+) \, [m^{3}s^{-1}]}$', vminlimit = 1e1)

    plot_compare(ts0, z0, nO2p/nNOp, n_ic2[:, 3, 1:]/n_ic2[:, 8, 1:], 'O2+/NO+ Density Ratio', r'$\mathrm{n(O_2^+)/n(NO^+)} \, [1]$')

    #no difference in major species
    # plot_compare(ts0, z0, nN2, n_ic2[:, 5, 1:], 'N2 Density', r'$\mathrm{n(N_2)} \, [m^{3}s^{-1}$]')
    # plot_compare(ts0, z0, nO2, n_ic2[:, 2, 1:], 'O2 Density', r'$\mathrm{n(O_2)} \, [m^{3}s^{-1}$]')
    # plot_compare(ts0, z0, nO , n_ic2[:, 1, 1:], 'O Density', r'$\mathrm{n(O)} \, [m^{3}s^{-1}$]')







# for i in range(15):
#     plot_compare(ts1[1:], z1, n_ic1[:, i, 1:], n_ic2[:, i, 1:], spstr[i], 'b', vminlimit=1)
# plt.show()


plt.show()
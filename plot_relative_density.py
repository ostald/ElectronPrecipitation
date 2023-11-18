import mat73
import matplotlib.pyplot as plt
import pickle
import numpy as np
import matplotlib.colors as colors

c = mat73.loadmat("/Users/ost051/Documents/PhD/Results/Comparison Flipchem-IonChem/flipchem_res_smooth.mat")
fc_iri = c['modelout']
[Tn_, Ti_, Te_, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = fc_iri.swapaxes(0,1)

c = mat73.loadmat("/Users/ost051/Documents/PhD/Results/Comparison Flipchem-IonChem/ElSpec_res.mat")["ElSpecQT_iqtOutliers_L5"]


f = "/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/2023.04.27_17_54_27_mixf=0/IC_res_14.pickle"
with open(f, 'rb') as pf:
     data = pickle.load(pf)

[ts_int, z_model, n_ic, eff_rr, ne_init, res] = data
n_ic = np.array([res[h].sol(c["te"]-c["te"][0]) for h in range(len(z_model))])
[e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(1, 0)

def plot_rel_abs_den(t, h, ne, ni, label):
    fig, ax = plt.subplots(nrows = 3, sharex = True)
    pc0 = ax[0].pcolormesh(t, h, ni/ne, vmin = 0)
    plt.colorbar(pc0, ax = ax[0], label = label + ' Rel. Density [1]')
    pc1 = ax[1].pcolormesh(t, h, ne)#, norm=colors.LogNorm())
    plt.colorbar(pc1, ax = ax[1], label = 'e Density [m-3]')
    pc2 = ax[2].pcolormesh(t, h, ni)#, norm=colors.LogNorm())
    plt.colorbar(pc2, ax = ax[2], label = label + ' Density [m-3]')
    fig.supxlabel('Time [s]')
    fig.supylabel('Height [km]')

tt = c["te"] - c["te"][0]
ii = 0

plot_rel_abs_den(tt, z_model, e[:, ii:], Op_4S[:, ii:], 'O+')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/Op_IC_RelAbundance.png')
plot_rel_abs_den(tt, z_model, c["neEnd"], nOp, 'O+')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/Op_ELSPEC_RelAbundance.png')


plot_rel_abs_den(tt, z_model, e[:, ii:], O2p[:, ii:], 'O2+')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/O2p_IC_RelAbundance.png')
plot_rel_abs_den(tt, z_model, c["neEnd"], nO2p, 'O2+')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/O2p_ELSPEC_RelAbundance.png')


plot_rel_abs_den(tt, z_model, e[:, ii:], NOp[:, ii:], 'NO+')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/NOp_IC_RelAbundance.png')
plot_rel_abs_den(tt, z_model, c["neEnd"], nNOp, 'NO+')
plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/NOp_ELSPEC_RelAbundance.png')
plt.show()
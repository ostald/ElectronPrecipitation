import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat

plt.rcParams.update({'font.size': 12})


eval_list = [
'2023.04.28_09_33_05_mixf=1',
'2023.05.03_14_07_04_mixf=0',
'2023.05.08_18_31_02_mixf=0',
'2023.05.10_17_18_42_mixf=0',
'2023.05.12_09_56_23_mixf=0',
'2023.04.27_17_54_27_mixf=0'
]

label_list = [
    'damped',
    'NO+/O2+ reversed',
    'NO+/O2+ = 1',
    'NO+ = O2+ = O+',
    '6*NO+ = 2*O2+ = O+',
    'IRI'
]

linestyles = [
    'dashed',
    'dashdot',
    #(5, (10, 3)),
    'dotted',
    (0, (3, 5, 1, 5)),
    (0, (5, 1)),
    'solid'
]

fig, axs = plt.subplots(nrows = 2, sharex=True, figsize=(6.4, 5))

for k, dir in enumerate(eval_list):
    f = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/'+ dir +'/IC_res_7.pickle'
    # mat = '/Users/ost051/Documents/PhD/ElectronPrecipitation/log/testing/'+ dir +'/ElSpec-iqt_IC_7.mat'
    #
    # spstr = ['e','O','O2','O2p','Np','N2','N2p','NO','NOp','H','Hp','O_1D','O_1S','N_2D','N_4S','O2p_a4P','Op_2D', \
    #          'Op_4S','Op_2P']
    #
    #
    # con = loadmat.loadmat(mat)["ElSpecOut"]
    # n_model = con["iri"]
    # [Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
    # eff_rr = con['alpha']
    # ts = con['ts']
    # ts = ts - ts[0]
    # z = con["h"]
    # ne = con["ne"]
    # q = con["q"]
    # Ie = con["Ie"]
    # E = con["E"][:-1]
    # [nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * ne
    # [ne_, Ti, Te, _] = con["par"].swapaxes(0, 1)

    with open(f, 'rb') as pf:
        data = pickle.load(pf)
        [ts, z, n_ic, eff_rr] = data[:4]

    i1 = 50
    if k==5 :
        axs[0].plot(ts[1:], n_ic[i1, 8, 1:]/n_ic[i1, 0, 1:], color = 'green' , linestyle = linestyles[k], label = 'NO+')
        axs[0].plot(ts[1:], n_ic[i1, 3, 1:]/n_ic[i1, 0, 1:], color = 'blue'  , linestyle = linestyles[k], label = 'O2+')
        axs[0].plot(ts[1:], n_ic[i1,17, 1:]/n_ic[i1, 0, 1:], color = 'orange', linestyle = linestyles[k], label = 'O+')
    else:
        axs[0].plot(ts[1:], n_ic[i1, 8, 1:] / n_ic[i1, 0, 1:], color='green' , linestyle=linestyles[k])
        axs[0].plot(ts[1:], n_ic[i1, 3, 1:] / n_ic[i1, 0, 1:], color='blue'  , linestyle=linestyles[k])
        axs[0].plot(ts[1:], n_ic[i1,17, 1:] / n_ic[i1, 0, 1:], color='orange', linestyle=linestyles[k])
    axs[0].text(0.1, 0.37, str(np.round(z[i1], 0)) + ' km', horizontalalignment='center', verticalalignment='center', transform=axs[0].transAxes)
    axs[0].plot(0, 0, 'x', alpha = 0)


    i2 = 18
    axs[1].plot([0, 0], [0, 0], color='black', linestyle=linestyles[k], label = str(k+1))# label = label_list[k])
    axs[1].plot(ts[1:], n_ic[i2, 8, 1:] / n_ic[i2, 0, 1:], color='green', linestyle=linestyles[k])
    axs[1].plot(ts[1:], n_ic[i2, 3, 1:] / n_ic[i2, 0, 1:], color='blue', linestyle=linestyles[k])
    axs[1].plot(ts[1:], n_ic[i2,17, 1:] / n_ic[i2, 0, 1:], color='orange', linestyle=linestyles[k])
    axs[1].text(0.1, 0.5, str(np.round(z[i2], 0)) + ' km', horizontalalignment='center', verticalalignment='center', transform=axs[1].transAxes)
    axs[1].plot(0, 0, 'x', alpha = 0)

    axs[1].set_xlabel('Time [s]')


fig.supylabel('Abundance [1]')

axs[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', mode="expand", borderaxespad=0., ncol = 3)
axs[1].legend(bbox_to_anchor=(0., 2.37, 1., .102), loc='lower left', mode="expand", borderaxespad=0., ncol = 6)
#axs[1].legend(bbox_to_anchor=(1.02, 0, 1.5, 1), loc='lower left', mode="expand", borderaxespad=0.)
#axs[1].legend()

plt.savefig('/Users/ost051/Documents/PhD/ElectronPrecipitation/writing/plots/convergence_bjorn.png')

plt.show()

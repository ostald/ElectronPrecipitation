import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat
import ionChem

direc = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.04.12_13_12_43 mixf=1/'
files = glob.glob(direc + '*.pickle')
print(sorted(files))

elspec_0 = loadmat.loadmat(direc + 'ElSpec-iqt_IC_0.mat')["ElSpecOut"]

matfiles = glob.glob(direc + 'ElSpec*.mat')

z = elspec_0["h"]
h = 0
species = 10
#spstr = ['e','O','Op','O2','O2p','N','Np','N2','N2p','NO','NOp','H','Hp']
spstr = ['e','O','O2','O2p','Np','N2','N2p','NO','NOp','H','Hp','O_1D','O_1S','N_2D','N_4S','O2p_a4P','Op_2D', \
         'Op_4S','Op_2P']


nax = 4

all_data = []

fig, axs = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    print(f)
    with open(f, 'rb') as pf:
        data = pickle.load(pf)
        all_data = [*all_data, data[:4]]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        #ne_init_old = data[4]

    f = direc + 'ElSpec-iqt_IC_' + str(i) + '.mat'
    elspec = loadmat.loadmat(f)["ElSpecOut"]
    ne_es = elspec["ne"]
    e_prod = elspec["q"]

    [e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, \
     Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(0, 1)

    axs.flat[i].plot(ts, e[0], label=r'$n_{e}$')
    axs.flat[i].plot(ts, ne_es[h, :], label=r'ElSpec $n_e$')
    axs.flat[i].plot(ts, np.sqrt(e_prod[0, :] / eff_rr[0, :]), '--', color='red', label=r'$ne_{ss}$')
    plt.legend(loc=2)
    plt.yscale('log')
    plt.xlabel('Time [s]')
    plt.ylabel(r'Density [m$^{-3}$]')
    # plt.title(c.name + ' Density')
    ax2 = axs.flat[i].twinx()
    ax2.plot(ts, e_prod[h, :], '.', color='green', label=r'$q_e$')
    ax2.set_yscale('log')
    ax2.legend(loc=1)
    ax2.set_ylabel(r'Electron production [m$^{-3}$s$^{-1}$]')
    # for t in ts_: plt.axvline(t, alpha = 0.1)
    #plt.tight_layout()
    #plt.show()
    #plt.savefig(direc + 'plots/' + c.name + '_Density_IC_' + str(iteration) + '.svg')
    # plt.savefig(direc + 'plots/' + c.name + '_Density_IC_' + str(iteration) + '.eps')

    axs.flat[i].set_title('Iteration ' + str(i))
fig.supxlabel('Time [s]')
fig.supylabel('Ratio of Charged Species')
axs.flat[0].legend(loc = 2)
fig.suptitle('Charged Species Stackplot at height index ' + str(h))


fig, axs = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    print(f)
    if True:
        data = all_data[i]
    #with open(f, 'rb') as pf:
        #data = pickle.load(pf)
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        #[ts, n_ic] = pickle.load(pf)
    [e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, \
     Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(0, 1) / n_ic.swapaxes(0, 1)[0]
    axs.flat[i].stackplot(ts, O2p[h], Np[h], N2p[h], NOp[h], Hp[h], O2p_a4P[h], Op_2D[h], Op_4S[h], Op_2P[h],
                 labels=['O2+', 'N+', 'N2+', 'NO+', 'H+', 'O2p_a4P', 'Op_2D', 'Op_4S', 'Op_2P'])
    axs.flat[i].set_title('Iteration ' + str(i))
fig.supxlabel('Time [s]')
fig.supylabel('Ratio of Charged Species')
axs.flat[0].legend(loc = 2)
fig.suptitle('Charged Species Stackplot at height index ' + str(h))

fig4, axs4 = plt.subplots()
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    if True:
        data = all_data[i]
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        #[ts, n_ic] = pickle.load(pf)
    [e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, \
     Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(0, 1)
    line = plt.plot(ts, n_ic.swapaxes(0, 1)[species][h], alpha=0.05, color='black', figure = fig4)
    plt.text(ts[-1], n_ic.swapaxes(0, 1)[species][h, -1], 'it'+str(i), color=line[0].get_color())
plt.yscale('log')
fig4.supxlabel('Time [s]')
fig4.supylabel('Number Density of ' + spstr[species])
fig4.suptitle('Number Density of ' + spstr[species] +' at height ' + str(z[h]))


fig31, axs31 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
vmin = 0
vmax = 0
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        den = n_ic.swapaxes(0, 1)[species]
        ma = np.max(den)
        if ma > vmax: vmax = ma
        mi = np.min(den)
        if mi > vmin: vmin = mi
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        #[ts, n_ic] = pickle.load(pf)
    den = n_ic.swapaxes(0, 1)[species]
    pc = axs31.flat[i].pcolormesh(ts, z, den, norm= mpl.colors.LogNorm(vmax=vmax))
    axs31.flat[i].set_title('Iteration ' + str(i))

fig31.suptitle('Number Density of ' + spstr[species])
fig31.supylabel('Number Density of ' + spstr[species])
fig31.supxlabel('Time [s]')
plt.colorbar(pc, ax = axs31)
#plt.show()

fig2, axs2 = plt.subplots()
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        #[ts, n_ic] = pickle.load(pf)
    [e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, \
     Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(0, 1)
    ratioO2pNOp = O2p/NOp
    line = plt.plot(ts, ratioO2pNOp[h], alpha=0.1, color='black', figure = fig2)
    plt.text(ts[-1], ratioO2pNOp[h, -1], 'it'+str(i), color=line[0].get_color())
fig2.supxlabel('Time [s]')
fig2.supylabel('Ratio O2+/NO+')
fig2.suptitle('Ratio O2+/NO+ at height index ' + str(h))


fig3, axs3 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
vmin = 0
vmax = 0
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        [e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, \
         Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(0, 1)
        ratioO2pNOp = O2p / NOp
        ma = np.max(ratioO2pNOp)
        if ma > vmax: vmax = ma
        mi = np.min(ratioO2pNOp)
        if mi > vmin: vmin = mi

for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        #[ts, n_ic] = pickle.load(pf)
    [e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, \
     Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(0, 1)
    ratioO2pNOp = O2p/NOp
    pc = axs3.flat[i].pcolormesh(ts, z, ratioO2pNOp, norm= mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    axs3.flat[i].set_title('Iteration ' + str(i))

fig3.suptitle('Ratio O2+/NO+')
fig3.supylabel('Altitude [km]')
fig3.supxlabel('Time [s]')
plt.colorbar(pc, ax = axs3)



fig4, axs4 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
pc4 = np.empty(len(files), dtype = 'object')
fig,ax = plt.subplots()

for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
    if i>0:
        # with open(direc + 'IC_res_'+str(i-1)+'.pickle', 'rb') as pf:
        #     data_o = pickle.load(pf)
        if True:
            data = all_data[i-1]
            ts_o = data[0][1:]
            z_o = data[1]
            n_ic_o = data[2][:, :, 1:]
            eff_rr_o = data[3][:, 1:]

    else: eff_rr_o = elspec_0["alpha"]
    d_effrr = eff_rr_o - eff_rr
    axs4.flat[i].set_title('Iteration ' + str(i))
    pc4 = axs4.flat[i].pcolor(ts, z, d_effrr, norm=mpl.colors.CenteredNorm(),
                    label='alpha', cmap='RdBu')
    fig4.colorbar(pc4, ax=axs4.flat[i])
    ax.plot(i, np.sum(np.abs(d_effrr.flat)), 'x', color = 'black')
    print(np.sum(np.abs(d_effrr.flat)))


fig4.suptitle('Deviation from previous eff. rec. rate [m3s-1]')
fig4.supylabel('Altitude [km]')
fig4.supxlabel('Time [s]')
#plt.colorbar(pc4, ax=axs4)
ax.set_ylabel(r"Sum of Deviations in $\alpha_{eff}$")
ax.set_xlabel("Iteration")
ax.set_yscale('log')

#plt.show()


fig41, axs41 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
pc41 = np.empty(len(files), dtype = 'object')
fig,ax = plt.subplots()

for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
    if i>0:
        if True:
            data = all_data[i-1]
            ts_o = data[0][1:]
            z_o = data[1]
            n_ic_o = data[2][:, :, 1:]
            eff_rr_o = data[3][:, 1:]
        #with open(direc + 'IC_res_'+str(i-1)+'.pickle', 'rb') as pf:
        #    [ts_o, z_o, n_ic_o, eff_rr_o] = pickle.load(pf)
    else: eff_rr_o = elspec_0["alpha"]
    d_effrr_r = (eff_rr_o - eff_rr) / eff_rr
    axs41.flat[i].set_title('Iteration ' + str(i))
    pc41 = axs41.flat[i].pcolor(ts, z, d_effrr_r, norm=mpl.colors.CenteredNorm(),
                    label='alpha', cmap='RdBu')
    fig41.colorbar(pc41, ax=axs41.flat[i])
    ax.plot(i, np.sum(np.abs(d_effrr_r.flat))/len(d_effrr_r.flat), 'x', color = 'black', label = 'Mean')
    ax.plot(i, np.max(np.abs(d_effrr_r.flat)), 'x', color = 'blue', label = 'Max')
    if i == 0: ax.legend()

fig41.suptitle('Relative Deviation from previous eff. rec. rate [m3s-1]')
fig41.supylabel('Altitude [km]')
fig41.supxlabel('Time [s]')
#plt.colorbar(pc4, ax=axs4)
ax.set_ylabel(r"Deviation in $\alpha_{eff}$")
ax.set_xlabel("Iteration")
ax.set_yscale('log')
#plt.show()
#exit()


fig5, axs5 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
vmin = 0
vmax = 0
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
        ma = np.max(eff_rr)
        if ma > vmax: vmax = ma
        mi = np.min(eff_rr)
        if mi > vmin: vmin = mi

for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    # with open(f, 'rb') as pf:
    #     data = pickle.load(pf)
    if True:
        data = all_data[i]
        ts = data[0][1:]
        z = data[1]
        n_ic = data[2][:, :, 1:]
        eff_rr = data[3][:, 1:]
    axs5.flat[i].set_title('Iteration ' + str(i))
    pc5 = axs5.flat[i].pcolor(ts, z, eff_rr, norm= mpl.colors.Normalize(vmin=vmin, vmax=vmax))

fig5.suptitle('Eff. Rec. Rate [m3s-1]')
fig5.supylabel('Altitude [km]')
fig5.supxlabel('Time [s]')
plt.colorbar(pc5, ax=axs5)

fig6, axs6 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
vmin = 0
vmax = 0
all_data_mat = []
for i, f in enumerate(matfiles[:nax**2]):
    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    print(f)
    elspec = loadmat.loadmat(f)["ElSpecOut"]
    all_data_mat = [*all_data_mat, elspec]
    ts = elspec["ts"]
    ne = elspec["ne"]
    ma = np.max(ne)
    if ma > vmax: vmax = ma
    mi = np.min(ne)
    if mi > vmin: vmin = mi

for i, f in enumerate(matfiles[:nax**2]):
    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    #elspec = loadmat.loadmat(f)["ElSpecOut"]
    elspec = all_data_mat[i]
    ts = elspec["ts"]
    ne = elspec["ne"]
    z = elspec["h"]
    axs6.flat[i].set_title('Iteration ' + str(i))
    pc6 = axs6.flat[i].pcolor(ts, z, ne, norm= mpl.colors.LogNorm(vmin=vmin, vmax=vmax))

fig6.suptitle('Electron Density Elspec')
fig6.supylabel('Altitude [km]')
fig6.supxlabel('Time [s]')
plt.colorbar(pc6, ax=axs6)


fig7, axs7 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
vmin = 0
vmax = 0
matfiles = glob.glob(direc + 'ElSpec*.mat')
for i, f in enumerate(matfiles[:nax**2]):
    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    #elspec = loadmat.loadmat(f)["ElSpecOut"]
    elspec = all_data_mat[i]
    ts = elspec["ts"]
    Ie = elspec["Ie"]
    ma = np.max(Ie)
    if ma > vmax: vmax = ma
    mi = np.min(Ie)
    if mi > vmin: vmin = mi

for i, f in enumerate(matfiles[:nax**2]):
    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    #elspec = loadmat.loadmat(f)["ElSpecOut"]
    elspec = all_data_mat[i]
    ts = elspec["ts"]
    Ie = elspec["Ie"]
    egrid = elspec["egrid"]
    axs7.flat[i].set_title('Iteration ' + str(i))
    pc7 = axs7.flat[i].pcolormesh(ts, egrid[:-1]/1e3, Ie, norm= mpl.colors.LogNorm(vmin = 1e5, vmax=vmax))
    axs7.flat[i].set_yscale('log')
    axs7.flat[i].set_ylim(1, 100)

fig7.suptitle('Electron Energy Spectrum')
fig7.supylabel('Energy [keV ]')
fig7.supxlabel('Time [s]')
plt.colorbar(pc7, ax=axs7)


fig8, axs8 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
vmin = 0
vmax = 0
matfiles = glob.glob(direc + 'ElSpec*.mat')
for i, f in enumerate(matfiles[:nax**2]):
    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    #elspec = loadmat.loadmat(f)["ElSpecOut"]
    elspec = all_data_mat[i]
    ts = elspec["ts"]
    Ie = elspec["Ie"]
    ma = np.max(Ie)
    if ma > vmax: vmax = ma
    mi = np.min(Ie)
    if mi > vmin: vmin = mi

for i, f in enumerate(matfiles[:nax**2]):
    f0 = direc + 'ElSpec-iqt_IC_0.mat'
    #elspec0 = loadmat.loadmat(f0)["ElSpecOut"]
    elspec0 = all_data_mat[0]
    Ie0 = elspec0["Ie"]

    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    # elspec = loadmat.loadmat(f)["ElSpecOut"]
    elspec = all_data_mat[i]
    ts = elspec["ts"]
    Ie = elspec["Ie"]
    egrid = elspec["egrid"]
    axs8.flat[i].set_title('Iteration ' + str(i))
    pc8 = axs8.flat[i].pcolormesh(ts, egrid[:-1]/1e3, ((Ie0-Ie)))#, norm= mpl.colors.LogNorm())#vmin = 1, vmax = 10e10))
    axs8.flat[i].set_yscale('log')
    axs8.flat[i].set_ylim(1, 100)

fig8.suptitle(' Difference in Electron Energy Spectrum to constant Ionospheric composotion')
fig8.supylabel('Energy [keV]')
fig8.supxlabel('Time [s]')
plt.colorbar(pc8, ax=axs8)
plt.show()

"""
lim = 1e5

f0 = sorted(matfiles)[0]
print('It 0: ', f0)
elspec0 = loadmat.loadmat(f0)["ElSpecOut"]
Ie0 = elspec0["Ie"]
Ie0_ = np.where(Ie0 < lim, np.ones(Ie0.shape)*lim, Ie0)

f1= sorted(matfiles)[-1]
print('It-1: ', f1)
elspec = loadmat.loadmat(f1)["ElSpecOut"]
Ie = elspec["Ie"]
Ie_ = np.where(Ie < lim, np.ones(Ie.shape)*lim, Ie)

ts = elspec["ts"]
egrid = elspec["egrid"]

fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
pc0 = axs.flat[0].pcolormesh(ts, egrid[:-1]/1e3, Ie0, norm=mpl.colors.LogNorm(vmin = 1e5))
pc1 = axs.flat[1].pcolormesh(ts, egrid[:-1]/1e3, Ie , norm=mpl.colors.LogNorm(vmin = 1e5))
# pc2 = axs.flat[2].pcolormesh(ts, egrid[:-1]/1e3, np.abs(Ie-Ie0), norm=mpl.colors.LogNorm(vmin = 1e5))#, label='alpha', cmap='RdBu')
# pc3 = axs.flat[3].pcolormesh(ts, egrid[:-1]/1e3, np.abs(Ie-Ie0)/Ie0_, norm=mpl.colors.LogNorm(vmin = 1, vmax = 1e3))#, label='alpha', cmap='RdBu')
# pc4 = axs.flat[4].pcolormesh(ts, egrid[:-1]/1e3, (Ie-Ie0)/Ie0_, norm=mpl.colors.CenteredNorm(halfrange = 10),
#                     label='alpha', cmap='RdBu')
# pc5 = axs.flat[5].pcolormesh(ts, egrid[:-1]/1e3, np.abs(Ie-Ie0)/(Ie0 + Ie))#, norm=mpl.colors.LogNorm(vmin = 1, vmax = 1e10))#, label='alpha', cmap='RdBu')
# pc6 = axs.flat[6].pcolormesh(ts, egrid[:-1]/1e3, Ie/Ie0_, norm=mpl.colors.LogNorm(vmin = 1, vmax = 1e3))#, label='alpha', cmap='RdBu')
# pc7 = axs.flat[7].pcolormesh(ts, egrid[:-1]/1e3, Ie_/Ie0, norm=mpl.colors.LogNorm(vmin = 1e-3, vmax = 1))#, label='alpha', cmap='RdBu')
pc8 = axs.flat[2].pcolormesh(ts, egrid[:-1]/1e3, Ie_/Ie0_, norm=mpl.colors.LogNorm(vmin = 1e-3, vmax = 1e3), cmap='RdBu')

fig.suptitle('Difference in Energy Spectrum')
fig.supylabel('Energy [km]')
fig.supxlabel('Time [s]')
plt.colorbar(pc0, ax = axs.flat[0])
plt.colorbar(pc1, ax = axs.flat[1])
plt.colorbar(pc8, ax = axs.flat[2])
#plt.colorbar(pc3, ax = axs.flat[3])
# plt.colorbar(pc4, ax = axs.flat[4])
# plt.colorbar(pc5, ax = axs.flat[5])
# plt.colorbar(pc6, ax = axs.flat[6])
# plt.colorbar(pc7, ax = axs.flat[7])
# plt.colorbar(pc8, ax = axs.flat[8])

axs.flat[0].set_yscale('log')
axs.flat[1].set_yscale('log')
axs.flat[2].set_yscale('log')
# axs.flat[3].set_yscale('log')
# axs.flat[4].set_yscale('log')
# axs.flat[5].set_yscale('log')
# axs.flat[6].set_yscale('log')
# axs.flat[7].set_yscale('log')
# axs.flat[8].set_yscale('log')

axs.flat[0].set_ylim(1, 100)
axs.flat[1].set_ylim(1, 100)
axs.flat[2].set_ylim(1, 100)
# axs.flat[3].set_ylim(1, 100)
# axs.flat[4].set_ylim(1, 100)
# axs.flat[5].set_ylim(1, 100)
# axs.flat[6].set_ylim(1, 100)
# axs.flat[7].set_ylim(1, 100)
# axs.flat[8].set_ylim(1, 100)
plt.show()



"""




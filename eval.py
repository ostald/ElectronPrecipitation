import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat
import ionChem

direc = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.02.16_17_58_37 mixf=0/'
files = glob.glob(direc + '*.pickle')

elspec_0 = loadmat.loadmat('/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/' +
                           'ElSpec-iqt_IC_0.mat')["ElSpecOut"]
z = elspec_0["h"]
h = 20
species = 10
#spstr = ['e','O','Op','O2','O2p','N','Np','N2','N2p','NO','NOp','H','Hp']
spstr = ['e','O','O2','O2p','Np','N2','N2p','NO','NOp','H','Hp','O_1D','O_1S','N_2D','N_4S','O2p_a4P','Op_2D', \
         'Op_4S','Op_2P']


nax = 3
fig, axs = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    print(f)
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
        #[ts, n_ic] = pickle.load(pf)
    [e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, \
     Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(0, 1) / n_ic.swapaxes(0, 1)[0]
    axs.flat[i].stackplot(ts, O2p[h], Np[h], N2p[h], NOp[h], Hp[h], O2p_a4P[h], Op_2D[h], Op_4S[h], Op_2P[h],
                 labels=['O+', 'O2+', 'N+', 'N2+', 'NO+', 'H+', 'O2p_a4P', 'Op_2D', 'Op_4S', 'Op_2P'])
    axs.flat[i].set_title('Iteration ' + str(i))
fig.supxlabel('Time [s]')
fig.supylabel('Ratio of Charged Species')
axs.flat[0].legend(loc = 2)
fig.suptitle('Charged Species Stackplot at height index ' + str(h))

fig4, axs4 = plt.subplots()
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
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
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
        den = n_ic.swapaxes(0, 1)[species]
        ma = np.max(den)
        if ma > vmax: vmax = ma
        mi = np.min(den)
        if mi > vmin: vmin = mi
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
        #[ts, n_ic] = pickle.load(pf)
    den = n_ic.swapaxes(0, 1)[species]
    pc = axs31.flat[i].pcolormesh(ts, z, den, norm= mpl.colors.LogNorm(vmax=vmax))
    axs31.flat[i].set_title('Iteration ' + str(i))

fig31.suptitle('Number Density of ' + spstr[species])
fig31.supylabel('Number Density of ' + spstr[species])
fig31.supxlabel('Time [s]')
plt.colorbar(pc, ax = axs31)


fig2, axs2 = plt.subplots()
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
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
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
        [e, O, O2, O2p, Np, N2, N2p, NO, NOp, H, Hp, O_1D, O_1S, N_2D, N_4S, O2p_a4P, \
         Op_2D, Op_4S, Op_2P] = n_ic.swapaxes(0, 1)
        ratioO2pNOp = O2p / NOp
        ma = np.max(ratioO2pNOp)
        if ma > vmax: vmax = ma
        mi = np.min(ratioO2pNOp)
        if mi > vmin: vmin = mi

for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
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
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
    if i>0:
        with open(direc + 'IC_res_'+str(i-1)+'.pickle', 'rb') as pf:
            [ts_o, z_o, n_ic_o, eff_rr_o] = pickle.load(pf)
    else: eff_rr_o = elspec_0["alpha"]
    d_effrr = eff_rr_o - eff_rr
    axs4.flat[i].set_title('Iteration ' + str(i))
    pc4 = axs4.flat[i].pcolor(ts, z, d_effrr, norm=mpl.colors.CenteredNorm(),
                    label='alpha', cmap='RdBu')
    fig4.colorbar(pc4, ax=axs4.flat[i])

fig4.suptitle('Deviation from previous eff. rec. rate [m3s-1]')
fig4.supylabel('Altitude [km]')
fig4.supxlabel('Time [s]')
#plt.colorbar(pc4, ax=axs4)


fig41, axs41 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
pc41 = np.empty(len(files), dtype = 'object')
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
    if i>0:
        with open(direc + 'IC_res_'+str(i-1)+'.pickle', 'rb') as pf:
            [ts_o, z_o, n_ic_o, eff_rr_o] = pickle.load(pf)
    else: eff_rr_o = elspec_0["alpha"]
    d_effrr_r = (eff_rr_o - eff_rr) / eff_rr
    axs41.flat[i].set_title('Iteration ' + str(i))
    pc41 = axs41.flat[i].pcolor(ts, z, d_effrr_r, norm=mpl.colors.CenteredNorm(),
                    label='alpha', cmap='RdBu')
    fig41.colorbar(pc41, ax=axs41.flat[i])

fig41.suptitle('Relative Deviation from previous eff. rec. rate [m3s-1]')
fig41.supylabel('Altitude [km]')
fig41.supxlabel('Time [s]')
#plt.colorbar(pc4, ax=axs4)



fig5, axs5 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
vmin = 0
vmax = 0
for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
        ma = np.max(eff_rr)
        if ma > vmax: vmax = ma
        mi = np.min(eff_rr)
        if mi > vmin: vmin = mi

for i, f in enumerate(files[:nax**2]):
    f = direc + 'IC_res_'+str(i)+'.pickle'
    with open(f, 'rb') as pf:
        [ts, z, n_ic, eff_rr] = pickle.load(pf)
    axs5.flat[i].set_title('Iteration ' + str(i))
    pc5 = axs5.flat[i].pcolor(ts, z, eff_rr, norm= mpl.colors.Normalize(vmin=vmin, vmax=vmax))

fig5.suptitle('Eff. Rec. Rate [m3s-1]')
fig5.supylabel('Altitude [km]')
fig5.supxlabel('Time [s]')
plt.colorbar(pc5, ax=axs5)

fig6, axs6 = plt.subplots(nrows=nax, ncols=nax, sharex=True, sharey=True)
vmin = 0
vmax = 0
matfiles = glob.glob(direc + 'ElSpec*.mat')
for i, f in enumerate(matfiles[:nax**2]):
    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    print(f)
    elspec = loadmat.loadmat(f)["ElSpecOut"]
    ts = elspec["ts"]
    ne = elspec["ne"]
    ma = np.max(ne)
    if ma > vmax: vmax = ma
    mi = np.min(ne)
    if mi > vmin: vmin = mi

for i, f in enumerate(matfiles[:nax**2]):
    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    elspec = loadmat.loadmat(f)["ElSpecOut"]
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
    elspec = loadmat.loadmat(f)["ElSpecOut"]
    ts = elspec["ts"]
    Ie = elspec["Ie"]
    ma = np.max(Ie)
    if ma > vmax: vmax = ma
    mi = np.min(Ie)
    if mi > vmin: vmin = mi

for i, f in enumerate(matfiles[:nax**2]):
    f = direc + 'ElSpec-iqt_IC_'+str(i)+'.mat'
    elspec = loadmat.loadmat(f)["ElSpecOut"]
    ts = elspec["ts"]
    Ie = elspec["Ie"]
    egrid = elspec["egrid"]
    axs7.flat[i].set_title('Iteration ' + str(i))
    pc7 = axs7.flat[i].pcolor(ts, egrid[:-1]/1e3, Ie, norm= mpl.colors.LogNorm(vmin = 1e5, vmax=vmax))
    axs7.flat[i].set_yscale('log')

fig7.suptitle('Electron Energy Spectrum')
fig7.supylabel('Altitude [km]')
fig7.supxlabel('Time [s]')
plt.colorbar(pc7, ax=axs7)
plt.show()


"""
plt.figure()
for r in model.all_reactions:
    #print(r.r_name)
    if r.r_name != 'gamma20 ': continue
    line = plt.plot(r.r_rate_t(Tn[:, 0], Ti[:, 0], Te[:, 0]), z_model)
    plt.text(r.r_rate_t(Tn[:, 0], Ti[:, 0], Te[:, 0])[0], z_model[0], r.r_name, color = line[0].get_color())
    
plt.xscale('log')
plt.title('Reaction Rate gamma20')
plt.ylabel('Altitude [km]')
plt.xlabel('Reaction Rate [m-3s-1]')
plt.show()
"""




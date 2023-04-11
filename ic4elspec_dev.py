import ic4elspec
import loadmat
import matplotlib.pyplot as plt
import numpy as np
import pickle

direc = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.04.05_18_42_52 mixf=1/'
file = 'ElSpec-iqt_IC_'
iteration = 0

chemistry_config = 'Data/other/Reaction rates full set.txt'

mixf = 0

#ic4elspec.ic(direc, chemistry_config, file, iteration, mixf = 0, test = True)

savedir = direc + "IC_res_" + str(iteration) + '.pickle'

with open(savedir, "rb") as f:
    data = pickle.load(f)
    eff_rr = data[3]
    res = data[5]

f = direc + file + str(iteration)
con = loadmat.loadmat(f)["ElSpecOut"]
h = con["h"]

plt.figure()
plt.plot(con["ne"][:, 0:20], h)
plt.figure()
plt.plot(con["alpha"][:, 0:42], h, label = 'Elsepc')
plt.plot(eff_rr[:, 0:42], h, label = 'IC')
plt.legend()
plt.figure()
plt.plot(con["q"][:, 0:20], h)
plt.plot(np.sqrt(con["q"]/con["alpha"])[:, 0:20], h)

plt.show()
quit()



print(eff_rr.shape)

plt.figure()
plt.plot(con["alpha"][:, 0], h, label = 'Elspec')
plt.plot(eff_rr[:,0], h, label = 'IC')
plt.legend()

n_ic = np.array([r.y for r in res])
o2p = n_ic[:, 3, :]
n2p = n_ic[:, 6, :]
nop = n_ic[:, 8, :]
ne_ic = n_ic[:, 0, :]

n_model = con["iri"]
ne = con["ne"]
[Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)

# normalise charged species to fit electron density
# necessary for every timestep, as ne changes in ElSpec, but n(ions) does not
[nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * ne

fig, ax = plt.subplots(1, 2)
ax.flat[0].plot(nO2p[:, 0]/ne[:, 0], h, label = 'Elspec O2p')
ax.flat[0].plot(o2p[:,0]/ne_ic[:, 0], h, label = 'IC')


ax.flat[1].plot(nNOp[:, 0]/ne[:, 0], h, label = 'Elspec NO')
ax.flat[1].plot(nop[:,0]/ne_ic[:, 0], h, label = 'IC')
ax.flat[1].legend()

# plt.figure()
# plt.plot(nOp[:, 0], h, label = 'Elspec')
# plt.plot(n2p[:,0], h, label = 'IC')
# plt.legend()


t = np.arange(-30*60, 0, 0.01)
o2p_t = res[0].sol(t)[3]/res[0].sol(t)[0]

ax.flat[0].plot(o2p_t[0], h[0], 'x', label = 't-30')
ax.flat[0].plot(o2p_t[-1], h[0], 'x', label = 't-0')
ax.flat[0].legend()

plt.figure()
plt.plot(t, o2p_t)
plt.show()

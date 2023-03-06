import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import loadmat
import ionChem
from matplotlib.animation import FuncAnimation

f = '/Users/ost051/Documents/PhD/Electron Precipitation/log/testing/2023.02.16_17_58_37 mixf=0/IC_res_8.pickle'

s = 3

spstr = ['e','O','O2','O2p','Np','N2','N2p','NO','NOp','H','Hp','O_1D','O_1S','N_2D','N_4S','O2p_a4P','Op_2D', \
         'Op_4S','Op_2P']

with open(f, 'rb') as pf:
    [ts, z, n_ic, eff_rr] = pickle.load(pf)

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = ax.plot([], [], 'ro')
txt = ax.text(0, n_ic.swapaxes(0, 1)[s].min(), 'h = ' + str(z[0]))


def init():
    ln, = plt.plot(ts, n_ic.swapaxes(0, 1)[s][0])
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Number Density')
    plt.yscale('log')
    ax.set_ylim([n_ic.swapaxes(0, 1)[s].min(), n_ic.swapaxes(0, 1)[s].max()])
    plt.title("Number density of " + spstr[s])
    return ln,

def update(frame):
    xdata = ts
    ydata = n_ic.swapaxes(0, 1)[s][frame]
    ln.set_data(xdata, ydata)
    txt.set_text('h = ' + str(z[frame]))
    return ln, txt,

ani = FuncAnimation(fig, update, frames=np.arange(len(z)),
                    init_func=init, blit=True)
plt.show()



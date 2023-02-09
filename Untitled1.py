#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import import_ipynb
import loadmat
import ionChemistry

ic4elspec.ic('/Users/ost051/Documents/PhD/ELSPEC-2022/Reactions Set Rees/',
             'ElSpec-iqt_IC_', 0)


# In[28]:


import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob

direc = '/Users/ost051/Documents/PhD/ELSPEC-2022/Reactions Set Rees/'
files = glob.glob(direc + '*.pickle')

h = 20

for f in files:
    with open(f, 'rb') as pf:
        res = pickle.load(pf)
    
    elspec_iri_ICsorted = np.array([r.y for r in res])
    [e,O,Op,O2,O2p,N,Np,N2,N2p,NO,NOp,H,Hp] = elspec_iri_ICsorted.swapaxes(0, 1)
    
    ts = res[0].t
    
    [re,rO,rOp,rO2,rO2p,rN,rNp,rN2,rN2p,rNO,rNOp,rH,rHp] = [e,O,Op,O2,O2p,N,Np,N2,N2p,NO,NOp,H,Hp]/e
    plt.figure()
    plt.stackplot(ts, rOp[h],rO2p[h],rNp[h],rN2p[h],rNOp[h],rHp[h],                   labels = ['O+', 'O2+', 'N+', 'N2+', 'NO+', 'H+'])
    plt.xlabel('Time [s]')
    plt.ylabel('Ratio of Charged Species')
    plt.legend(loc = 2)
    plt.title('Charged Species Stackplot at height index ' + str(h))


# In[32]:


elspec_0 = loadmat.loadmat('/Users/ost051/Documents/PhD/ELSPEC-2022/'+                           'Reactions Set Rees/ElSpec-iqt_IC_0.mat')["ElSpecOut"]
n_model = elspec_0["iri"]
n_model.shape
[Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp] = n_model.swapaxes(0, 1)


# In[32]:


import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import import_ipynb
import loadmat
#from labellines import labelLine, labelLines


elspec_0 = loadmat.loadmat('/Users/ost051/Documents/PhD/ELSPEC-2022/'+                    'Reactions Set Rees/ElSpec-iqt_IC_0.mat')["ElSpecOut"]
n_model = elspec_0["iri"]
[Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp] = n_model.swapaxes(0, 1)

direc = '/Users/ost051/Documents/PhD/ELSPEC-2022/Reactions Set Rees/'
files = glob.glob(direc + '*.pickle')


h = 0

fig = plt.figure()

for i, f in enumerate(files):
    with open(f, 'rb') as pf:
        res = pickle.load(pf)
    
    elspec_iri_ICsorted = np.array([r.y for r in res])
    [e,O,Op,O2,O2p,N,Np,N2,N2p,NO,NOp,H,Hp] = elspec_iri_ICsorted.swapaxes(0, 1)
    
    ts = res[0].t
    
    [re,rO,rOp,rO2,rO2p,rN,rNp,rN2,rN2p,rNO,rNOp,rH,rHp] = [e,O,Op,O2,O2p,N,Np,N2,N2p,NO,NOp,H,Hp]/e

    ratioO2pNOp = O2p/NOp
    
    #if i==0: plt.plot(ts, ratioO2pNOp[h], label = 'i0')
    #else: pass
    #if f==files[-1]: pass
    #else: 
    line = plt.plot(ts, ratioO2pNOp[h], alpha = 0.1, color = 'black')
    #line = plt.plot(lossRate, self.z_model/1e3, label = c.name)
    plt.text(ts[-1], ratioO2pNOp[h, -1], i, color = line[0].get_color())
    

    
line = plt.plot(ts, (nO2p/nNOp)[h], label = 'i0')   
plt.text(ts[-1], (nO2p/nNOp)[h, -1], 0, color = line[0].get_color())
plt.xlabel('Time [s]')
plt.ylabel('Ratio O2p/NO2')
plt.legend()
plt.title('Charged Species Stackplot at height index ' + str(h))
#plt.yscale('log')


fig2 = plt.figure()
f = files[3]
with open(f, 'rb') as pf:
    res = pickle.load(pf)
    
elspec_iri_ICsorted = np.array([r.y for r in res])
[e,O,Op,O2,O2p,N,Np,N2,N2p,NO,NOp,H,Hp] = elspec_iri_ICsorted.swapaxes(0, 1)
    
ts = res[0].t

[re,rO,rOp,rO2,rO2p,rN,rNp,rN2,rN2p,rNO,rNOp,rH,rHp] = [e,O,Op,O2,O2p,N,Np,N2,N2p,NO,NOp,H,Hp]/e
ratioO2pNOp = O2p/NOp

pc = plt.pcolormesh(ts, elspec_0["h"], ratioO2pNOp)
plt.colorbar(pc)


# In[27]:


import import_ipynb
import loadmat
import ionChem
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'widget')

elspec_0 = loadmat.loadmat('/Users/ost051/Documents/PhD/ELSPEC-2022/'+                    'Reactions Set Rees/ElSpec-iqt_IC_0.mat')["ElSpecOut"]
n_model = elspec_0["iri"]
[Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp] = n_model.swapaxes(0, 1)
ts = elspec_0["ts"]


chemistry_config = '/Users/ost051/Documents/PhD/ELSPEC-2022/'+                        'Reactions Set Rees/Reaction rates.txt'
z_model = elspec_0["h"]    
model = ionChem.ionChem(chemistry_config, z_model)

plt.figure()
plt.pcolormesh(ts, z_model, (Ti + Tn)/2)
plt.colorbar()

plt.figure()
for r in model.all_reactions:
    #print(r.r_name)
    if r.r_name != 'gamma20 ': continue
    line = plt.plot(r.r_rate_t(Tn[:, 0], Ti[:, 0], Te[:, 0]), z_model)
    plt.text(r.r_rate_t(Tn[:, 0], Ti[:, 0], Te[:, 0])[0], z_model[0], r.r_name, color = line[0].get_color())
    
plt.xscale('log')


# In[75]:


line[0].get_data()


# In[68]:


(4000/5778)**4 *25**2


# In[25]:


import numpy as np
Tr = np.arange(0, 2000, 100)
Ti = np.arange(0, 2000, 100)
Tn = np.arange(0, 2000, 100)
print(model.all_reactions[17].r_rate_string)
r = eval(model.all_reactions[17].r_rate_string)
print(r)
plt.figure()
plt.plot(Tr, r)
plt.yscale('log')


# In[ ]:





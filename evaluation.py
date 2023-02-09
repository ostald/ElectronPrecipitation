#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Loading stored results and processing


# In[11]:


# to dump results:
# setup_.datadump(model)

import pickle

def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break
                
data = list(loadall(log_directory + '/Simulation/log.p'))

data[0].all_species[3].density


# In[12]:


time = np.array([i.time for i in data])
densities = np.array([[j.density for j in i.all_species] for i in data])


# In[13]:


get_ipython().run_line_magic('matplotlib', 'widget')
plt.figure()
for c in data[0].ions:
    i = c.c_ID
    plt.plot(time, densities[:, i, -1], label=c.name)
plt.yscale('log')
plt.legend()


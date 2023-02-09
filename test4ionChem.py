#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import import_ipynb

import organizing
import loadMSIS
import loadIRI
import EnergyDissipation 
import loadEISCAT

'''
This is the main file. It is used to start the simulation.
Configuration files from IRI and MSIS models are specified here.
Eiscat data files are specified here.
'''

msis_config = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/msis.txt'
iri_config  = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/iri.txt'
chemistry_config = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/Reaction rates.txt'
path_eiscat_data = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Data/'

setup_ = organizing.setup(msis_config, iri_config, chemistry_config, path_eiscat_data)
log_directory = setup_._log_directory

[z_msis
 , n_o1_msis
 , n_n2_msis
 , n_o2_msis
 , mass_density
 , temp_n_msis
 , scale_height_msis] = loadMSIS.loadMSIS(log_directory + '/Config/msis.txt')

[z_iri     
 , ne_iri    
 , t_neutral 
 , t_ion     
 , t_e       
 , rel_o_p   
 , rel_n_p   
 , rel_h_p   
 , rel_he_p  
 , rel_o2_p  
 , rel_no_p] = loadIRI.loadIRI(log_directory + '/Config/iri.txt')


# In[2]:


log_directory


# In[3]:


[z_radar_time
 , parameters
 , time
 , az
 , el
 , r_errors ] = loadEISCAT.loadEISCATdata(path_eiscat_data, ['r_h', 'r_param', 'r_time', 'r_az', 'r_el', 'r_error'])

z_radar = z_radar_time[0] #for static measurements, the information in r_h stays the same for all measurement points.
ne_meas = parameters[:, :, 0]
dne_meas = r_errors[:, :, 0]
vel = parameters[:, :, 4]
vel_err = r_errors[:, :, 4]

#conversion to SI units:
z_radar = z_radar*1e3
r_time = time[:,0,3]*3600 + time[:,0,4]*60 + time[:,0,5]

del parameters, r_errors, z_radar_time, az, el, time


# In[4]:


#differential flux: energy binning: 
e_bin_lims = np.logspace(2, 2+3, 501)

#defining parameters for electron production model
inc_B = 0
z_min = 70000
z_max = 200000
res_model = 100
z_sd_min = 0
z_sd_max = 800e3
res_sd = 10

e_prod, z_model = EnergyDissipation.e_production(z_min
                                                 , z_max
                                                 , res_model
                                                 , z_sd_min
                                                 , z_sd_max
                                                 , res_sd
                                                 , e_bin_lims
                                                 , z_msis
                                                 , mass_density
                                                 , inc_B)


# In[6]:


#setup of ion chemistry
import ionChem

model = ionChem.ionChem('/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/Reaction rates test.txt'
, np.arange(10))

#model.assign_densities(z_model
#                     , z_msis
#                     , n_o1_msis
#                     , n_n2_msis
#                     , n_o2_msis
#                     , z_iri
#                     , ne_iri        
#                     , rel_o_p   
#                     , rel_n_p   
#                     , rel_h_p   
#                     , rel_he_p  
#                     , rel_o2_p  
#                     , rel_no_p)


z_model = np.arange(0, 10)
model.all_species[0].density  = np.ones(10)*1.1
model.all_species[1].density  = np.zeros(10)
model.all_species[2].density  = np.ones(10)*0.9
model.all_species[3].density  = np.ones(10)*2



# In[7]:


import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'widget')

def plot_snapshot(model):
    plt.figure()
    for c in model.all_species:
        plt.plot(c.density, z_model/1e3, label = c.name)
    #plt.xscale('log')
    plt.xlabel('Density [m-3]')
    plt.ylabel('Altitude [km]')
    plt.title('Density profiles')
    plt.legend()
    plt.xlim(-0.5,2.5)


# In[8]:


dt = 0.1
t = 0
setup_.datadump(model)


# In[9]:


for i in range(100):
    model.update_timestep(dt)
    t = t + dt
    #plot_snapshot(model)
    setup_.datadump(model)


# In[ ]:


import pickle

def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break
                
data = list(loadall(log_directory + '/Simulation/log.p'))

data[1].all_species[3].density


# In[11]:


time = np.array([i.time for i in data])
densities = np.array([[j.density for j in i.all_species] for i in data])


# In[12]:


plt.figure()
plt.plot(time, densities[:, 0, 0], 'x')
plt.plot(time, densities[:, 1, 0], 'x')
plt.plot(time, densities[:, 2, 0], 'x')
plt.plot(time, densities[:, 3, 0], 'x')


# In[16]:


type(model.all_species[0].charge)


# In[14]:


any(np.array([1, 2, 3])>2.9)


# In[20]:


[i for i in [1, 2, 3] if i>2]


# In[ ]:





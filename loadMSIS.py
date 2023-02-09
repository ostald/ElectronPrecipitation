#!/usr/bin/env python
# coding: utf-8

# In[33]:


import numpy as np

import import_ipynb
from organizing import pcfg
con = pcfg

_printing = con.print

_kb = 1.380649e-23 #in SI units 
_re = 6371000      #in m
_g = 9.81          #in m/s^2

def loadMSIS_new(file):
    #reading msis data
    #get length, define data array size
    f = open(file, 'r')
    n_lines = len(f.readlines())
    data = np.zeros([n_lines - 45, 9])
    f.close()
    
    #get data into data array
    f = open(file, 'r')
    for i, line in enumerate(f.readlines()[45:-1]):
        data[i] = line.replace('\\', '').split()
    f.close()
        
    #print first few lines
    if con.print:
        f = open(file, 'r')
        for i in range(25):
            print(f.readline().replace('\n', ''))
            continue
        f.close()
    
    #more intuitive names and SI unit conversion
    z_msis       = data[:, 0]*1e3      #[height in m]
    n_o1_msis    = data[:, 1]*1e6      #[number density in m^-3]
    n_n2_msis    = data[:, 2]*1e6      #[number density in m^-3]
    n_o2_msis    = data[:, 3]*1e6      #[number density in m^-3]
    mass_density = data[:, 4]*1e3      #[mass density (density in kg/m^3]
    temp_n_msis  = data[:, 5]          #[neutral temperature in K]
    n_h_msis     = data[:, 7]*1e6      #[number density in m^-3]
    n_n_msis     = data[:, 8]*1e6      #[number density in m^-3]
    
    number_density = np.sum(np.array([n_o1_msis, n_n2_msis, n_o2_msis]), 0)
    mean_molecular_mass = mass_density/number_density
    
    avarage_molecular_mass = mass_density/number_density
    
    #calculate scale height:
    global scale_height_msis
    scale_height_msis = _kb * temp_n_msis * (_re**2 + z_msis**2) / (avarage_molecular_mass * _g * _re**2)
    
    return [z_msis, n_o1_msis, n_n2_msis, n_o2_msis, n_h_msis, n_n_msis, mass_density, temp_n_msis, scale_height_msis]


# In[6]:


def loadMSIS(file):
    #reading msis data
    #get length, define data array size
    f = open(file, 'r')
    n_lines = len(f.readlines())
    data = np.zeros([n_lines - 15, 7])
    f.close()
    
    #get data into data array
    f = open(file, 'r')
    for i, line in enumerate(f.readlines()[15:]):
        data[i] = line.split()
    f.close()
        
    #print first few lines
    if con.print:
        f = open(file, 'r')
        for i in range(25):
            print(f.readline().replace('\n', ''))
            continue
        f.close()
    
    #more intuitive names and SI unit conversion
    z_msis       = data[:, 0]*1e3      #[height in m]
    n_o1_msis    = data[:, 1]*1e6      #[number density in m^-3]
    n_n2_msis    = data[:, 2]*1e6      #[number density in m^-3]
    n_o2_msis    = data[:, 3]*1e6      #[number density in m^-3]
    mass_density = data[:, 4]*1e3      #[mass density (density in kg/m^3]
    temp_n_msis  = data[:, 5]          #[neutral temperature in K]
    
    number_density = np.sum(np.array([n_o1_msis, n_n2_msis, n_o2_msis]), 0)
    mean_molecular_mass = mass_density/number_density
    
    avarage_molecular_mass = mass_density/number_density
    
    #calculate scale height:
    global scale_height_msis
    scale_height_msis = _kb * temp_n_msis * (_re**2 + z_msis**2) / (avarage_molecular_mass * _g * _re**2)
    
    return [z_msis, n_o1_msis, n_n2_msis, n_o2_msis, mass_density, temp_n_msis, scale_height_msis]


# In[21]:


#visualization
if con.print:
    import matplotlib.pyplot as plt
    get_ipython().run_line_magic('matplotlib', 'widget')
    from scipy.interpolate import CubicSpline, interp1d, PchipInterpolator

    MSISFile = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/msis.txt'
    [z_msis
     , n_o1_msis
     , n_n2_msis
     , n_o2_msis
     , mass_density
     , temp_n_msis
     , scale_height_msis] = loadMSIS(MSISFile)
    
    plt.figure()
    plt.plot(n_o1_msis, z_msis/1e3, label = 'O')
    plt.plot(n_o2_msis, z_msis/1e3, label = 'O2')
    plt.plot(n_n2_msis, z_msis/1e3, label = 'N2')
    plt.ylabel('Height [km]')
    plt.xlabel('Density [m-3]')
    plt.xscale('log')
    plt.legend()


# In[7]:


if con.print:
    plt.figure()
    plt.plot(mass_density, z_msis/1e3, 'x-', label = 'mass:dens')
    z = np.arange(70000, 200000, 1)
    int_lin  = np.exp(interp1d(z_msis[:2], np.log(mass_density[:2]), fill_value = 'extrapolate')(z))
    int_cs = np.exp(CubicSpline(z_msis, np.log(mass_density))(z))
    
    #plt.plot(int_cs, z/1e3)
    plt.plot(int_lin, z/1e3)
    plt.plot(np.array([*int_lin[z<80000], *int_cs[z>=80000]]), z/1e3)
    
    print(int_cs[z>=80000])
    
    plt.ylabel('Height [km]')
    plt.xlabel('Density [m-3]')
    plt.xscale('log')
    plt.legend()


# In[17]:


if con.print:
    plt.figure()
    plt.plot((mass_density), z_msis/1e3, 'x', label = 'mass:dens')
    z = np.arange(70000, 200000, 1)
    plt.plot(np.exp(interp1d(z_msis, np.log(mass_density), fill_value = 'extrapolate')(z)), z/1e3, label = 'exp(lin(log)))')
    plt.plot(np.exp(CubicSpline(z_msis, np.log(mass_density))(z)), z/1e3, label = 'exp(cs(log)))')
    plt.plot(np.exp(PchipInterpolator(z_msis, np.log(mass_density))(z)), z/1e3, label = 'exp(pchip(log))')
    plt.ylabel('Height [km]')
    plt.xlabel('Density [m-3]')
    plt.xscale('log')
    plt.legend()


# In[ ]:





# In[ ]:





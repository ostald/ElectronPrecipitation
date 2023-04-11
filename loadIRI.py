#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np

import import_ipynb
from organizing import pcfg
con = pcfg

printing = con.print
printing = 1

def loadIRI(iri_file):
    #iri_file = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/iri.txt'
    skip_lines = 27
    
    #get length, define data array size
    f = open(iri_file, 'r')
    n_lines = len(f.readlines())
    data = np.zeros([n_lines - skip_lines, 15])
    f.close()
    
    #get data into data array
    f = open(iri_file, 'r')
    for i, line in enumerate(f.readlines()[skip_lines:]):
        data[i] = line.split()
    f.close()
        
    #print first few lines
    if con.print:
        f = open(iri_file, 'r')
        for i in range(42):
            print(f.readline().replace('\n', ''))
            continue
        f.close()
        
    #more intuitive names and SI unit conversion
    z_iri        = data[:, 0]*1e3         #[height in m]
    ne_iri       = data[:, 1]*1e6         #[number density in m^3]
    t_neutral    = data[:, 3]             #[Neutral Temperature in K]
    t_ion        = data[:, 4]             #[Ion Temperature in K]
    t_e          = data[:, 5]             #[Electron Temperature in K]
    rel_o_p      = data[:, 6]/100             #[O+ percentage of ions]
    rel_n_p      = data[:, 7]/100             #[N+ percentage of ions]
    rel_h_p      = data[:, 8]/100             #[H+ percentage of ions]
    rel_he_p     = data[:, 9]/100             #[He+ percentage of ions]
    rel_o2_p     = data[:,10]/100             #[O2+ percentage of ions]
    rel_no_p     = data[:,11]/100             #[NO+ percentage of ions]
    
#    #replace all values where t_e is -1 to 300:
#    t_e       = np.array([t if t!= -1 else 300 for t in t_e])
#    t_ion     = np.array([t if t!= -1 else 300 for t in t_ion])
#    t_neutral = np.array([t if t!= -1 else 300 for t in t_neutral])
    
    return [z_iri     
            , ne_iri    
            , t_neutral 
            , t_ion     
            , t_e       
            , rel_o_p   
            , rel_n_p   
            , rel_h_p   
            , rel_he_p  
            , rel_o2_p  
            , rel_no_p
            ]


# In[2]:


if printing:
    iri_file = '/Users/ost051/Documents/PhD/Electron Precipitation/Data/other/iri.txt'

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
        , rel_no_p
     ] = loadIRI(iri_file)
    
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.plot(t_neutral, z_iri/1e3, label = 'Tn')
    plt.plot(t_ion, z_iri/1e3, label = 'Ti')
    plt.plot(t_e, z_iri/1e3, label = 'Te')
    plt.xlabel('Tempereature [K]')
    plt.ylabel('Altitude [km]')
    
    plt.figure()
    plt.plot(rel_o_p, z_iri/1e3, label = 'O+')
    plt.plot(rel_n_p, z_iri/1e3, label = 'N+')
    plt.plot(rel_h_p, z_iri/1e3, label = 'H+')
    plt.plot(rel_he_p, z_iri/1e3, label = 'He+')
    plt.plot(rel_o2_p, z_iri/1e3, label = 'O2+')
    plt.plot(rel_no_p, z_iri/1e3, label = 'NO+')
    plt.legend()
    plt.xlabel('Percent [%]')
    plt.ylabel('Altitude [km]')
    
    plt.figure()
    plt.plot(ne_iri, z_iri/1e3)
    plt.xlabel('Electron Density [m-3]')
    plt.ylabel('Altitude [km]')
    
    plt.figure()
    plt.plot(t_e, z_iri/1e3)
    plt.xlabel('ELectron Temperature [K]')
    plt.ylabel('Altitude [km]')
    plt.show()

    plt.figure()
    import import_ipynb
    from loadMSIS import *
    MSISFile = '/Users/ost051/Documents/PhD/Electron Precipitation/Data/other/msis.rtf'
    [z_msis, n_o1_msis, n_n2_msis, n_o2_msis, mass_density, temp_n_msis, scale_height_msis] = loadMSIS(MSISFile)
    plt.plot(n_o1_msis, z_msis/1e3, label = 'O')
    plt.plot(n_o2_msis, z_msis/1e3, label = 'O2')
    plt.plot(n_n2_msis, z_msis/1e3, label = 'N2')
    plt.plot(rel_o_p*ne_iri, z_iri/1e3, label = 'O+')
    plt.plot(rel_n_p*ne_iri, z_iri/1e3, label = 'N+')
    plt.plot(rel_h_p*ne_iri, z_iri/1e3, label = 'H+')
    plt.plot(rel_he_p*ne_iri, z_iri/1e3, label = 'He+')
    plt.plot(rel_o2_p*ne_iri, z_iri/1e3, label = 'O2+')
    plt.plot(rel_no_p*ne_iri, z_iri/1e3, label = 'NO+')
    plt.ylabel('Height [km]')
    plt.xlabel('Density [m-3]')
    plt.xscale('log')
    plt.legend()


# In[ ]:





# In[ ]:





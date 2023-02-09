#!/usr/bin/env python
# coding: utf-8

# In[1]:


#functions related to the radar


# In[ ]:





# In[9]:


def pulse_smearer(pulse_height, resolution, electron_profile):
    """
    Simulates smearing of observables (such as the electron profile) due to finite pulse length
    use a gliding average on the production profiles q (= fm[:, i])
    """
    pulse_height = 2000 #pulse height in m
    ps_range = pulse_height / res_model
    pulse_smearer = np.ones(np.int_(ps_range))
    pulse_smearer = pulse_smearer / np.sum(pulse_smearer) #ensures that the pulse smearer always has norm 1
    
    e_prod_ps = np.array([np.convolve(height_prof, pulse_smearer, 'same') for height_prof in electron_profile.T]).T
    return e_prod_ps


# In[ ]:


def ne_height_int(z_model, ne_tp1, z0, z1):
    """
    Calculates the mean electron density as the radar would measure it,
    from a high resolved height profile (z_fm, ne_tp1) 
    to one measurement between [z0, z1]
    
    Parameters:
    z_fm: list of heights along which the  electron density is known
        ndarray
    ne_tp1: list of electron density along the heights given by z_fm
        ndarray
    z0: lower integration boundary
        scalar
    z1: upper integration boundary
        scalar
        
    Returns:
    ne_int: The height integrated electron density, normalized by the difference z0-z1
    """
    #if (ne_tp1<0).any(): raise RuntimeError('Negative Electron Density')
    res = 1 # resolution in meters
    z = np.arange(z0, z1, res)
    ne_ip = interp1d(z_model, ne_tp1)(z)
    ne_int = (np.sum(ne_ip) * res + ne_ip[-1] * (z1 - (z[-1]+res))) / (z1 -z0)
    #if (ne_int<0).any(): raise RuntimeError('Negative Electron Density')
    return ne_int

ne_int = np.zeros(len(z_radar))
for i in range(len(z_radar)):
    ne_int[i] = ne_height_int(z_model, ne_iri_model, zr_lim[i], zr_lim[i+1])
    
if printing:
    plt.figure()
    plt.plot(ne_iri_fm, z_fm/1e3, label = 'high res')
    plt.plot(ne_int, z_r/1e3, 'x', label = 'radar res')
    plt.legend()


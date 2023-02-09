#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import import_ipynb
import warnings

from scipy.interpolate import CubicSpline, interp1d
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import AuroralPrecipitation as ap
from organizing import pcfg
con = pcfg


# In[2]:


# a collection of useful functions

def diff_flux_kev(e, p):
    """
    Calculates the differential flux as  E*exp(polynomial in E) for energies e
    with p the factors of a Polynomial of third degree
    
    Parameters:
    e: energy to evaluate the function [ev]
        scalar or ndarray
    p: coefficients as a list 
       IN [keV^-n] TO SCALE THE VALUES OF P CLOSER TO 1!
        ndarray
    
    Returns:
    diff_flux: the differential energy flux in [ev-1 m-2 s-1]
        ndarray
    """
    ke = e/1e3
    diff_flux = e * np.exp(p[0] + p[1]*ke + p[2]*ke**2 + p[3]*ke**3)
    return diff_flux


def pulse_smearing_4_e_prod(e_prod, pulse_height = 2000, res = 100):
    """
    Simulates the pulse smearing effect of a radar measurement by convolving the production profiles with a filter.
    The filter as of now is a box function.
    
    Parameters:
    e_prod: The matrix a of q = A \phi, i.e. the height profiles of the electron production rates due to a electron energy spectrum \phi.
        ndarray
    pulse_height: the height extend of a pulse
        scalar
    res: the height resolution of e_prod
        scalar
    
    Returns:
    e_prod_ps: the pulse-smaered electron production matrix
    """
    ps_range = pulse_height / res
    i_correction = np.int_(np.ceil(ps_range/2)) #to correct for teh last values, falling off to 0 after ps
    pulse_smearer = np.ones(np.int_(ps_range))
    pulse_smearer = pulse_smearer / np.sum(pulse_smearer) #ensures that the pulse smearer always has norm 1

    e_prod_ps = np.array([np.convolve(height_prof, pulse_smearer, 'same') for height_prof in e_prod.T]).T
    e_prod_ps[i_correction:, :] = e_prod[i_correction:, :] #to correct for the last values, falling off to 0 after ps
    return e_prod_ps
        
    
def ne_height_int(z, ne, z0, z1):
    """
    Calculates the mean electron density as the radar would measure it,
    from a high resolved height profile (z, ne) 
    to one measurement between [z0, z1]
    (needs to be called for every height interval seperately).
    
    Parameters:
    z: list of heights along which the  electron density is known
        ndarray
    ne: list of electron density along the heights given by z
        ndarray
    z0: lower integration boundary
        scalar
    z1: upper integration boundary
        scalar
        
    Returns:
    ne_int: The height integrated electron density, normalized by the difference z0-z1
    """
    if (ne<0).any():
        print(ne)
        raise RuntimeError('Negative Electron Density')
    res = 1 # resolution for interpolation in in meters
    z_new = np.arange(z0, z1, res)            
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ne_ip = np.exp(interp1d(z, np.log(ne), fill_value = 'extrapolate')(z_new)) #double check!!

    #ne_ip = np.exp(interp1d(z, np.log(ne), fill_value = 'extrapolate')(z_new)) #double check!!
    #if con.print == 1: plt.plot(ne, z/1e3, 'x', color = 'tab:blue')
    #if con.print == 1: plt.plot(ne_ip, z_new/1e3, color = 'tab:orange')
    ne_ip = interp1d(z, ne, fill_value = 'extrapolate')(z_new)
    ne_int = (np.sum(ne_ip) * res + ne_ip[-1] * (z1 - (z_new[-1]+res))) / (z1 -z0)
    if (ne_int<0).any():
        print('Error raised after height integration', z0, z1, ne_int)
        raise RuntimeError('Negative Electron Density')
    return ne_int



def find_ne0(param, e_bin_mean, ne_meas, dne_meas, e_prod_ps, eff_recombinationRate, z_model, zr_lim):
    """
    Used for calculating the first estimate of ne.
    Uses steady-state assumption i.e. ne = sqrt(q/eff_rr)
    """
    
    def sumofsq_ne0(par4diff_flux, e_bin_mean, ne_obs, dne_obs, e_prod_ps, eff_recombinationRate, z_model, zr_lim):
        diff_flux = ap.diff_flux_kev(e_bin_mean, par4diff_flux)
        ne_mod    = np.sqrt(np.dot(e_prod_ps, diff_flux)/eff_recombinationRate)
        ne_mod_r  = np.array([ap.ne_height_int(z_model, ne_mod, zr_lim[i], zr_lim[i+1]) for i in range(len(zr_lim)-1)])
        sumofsq   = np.sum((ne_obs - ne_mod_r)**2/dne_obs**2)
        return sumofsq
    
    #minimize
    from scipy.optimize import minimize
    result_t0 = minimize(sumofsq_ne0, 
                         param,
                         method = 'Nelder-Mead',
                         bounds = [(None, None),(None, None),(-np.inf, np.inf),(-np.inf, 0)], 
                         args=(e_bin_mean, ne_meas, dne_meas, e_prod_ps, eff_recombinationRate, z_model, zr_lim)
                         )
    return result_t0
 


    


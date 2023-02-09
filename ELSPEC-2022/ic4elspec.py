#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scipy.io as spio
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
import import_ipynb
import ionChem
import numpy as np
import matplotlib.pyplot as plt


def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def ic(file):
    mat = loadmat(file)
    con = mat["ElSpecOut"]

    n_model = con["iri"]
    [Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp] = n_model[:, :, 0].T
    temp = n_model[:, :3, :].T
    ts = con["ts"]
    te = con["te"] - ts[0]
    ts = ts - ts[0]
    e_prod = con["q"].T
    
    def stepped_prod_t(prod, t):
    """
    returns the production according to the ELSPEC model
    of any species where prod is defined
    at arbitrary times
    using the last ts from ELSPEC
    """
    if t<ts[0]:
        return 0
    if t >te[-1]:
        return 0
    else:
        i_max_ts = len(ts[ts<=t])-1
        prod_t = prod[i_max_ts]
        return prod_t
    
    chemistry_config = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/Reaction rates.txt'
    z_model = con["h"]
    
    model = ionChem.ionChem(chemistry_config, z_model)
    
    #assign densities
    model.N2.density  = nN2
    model.O2.density  = nO2
    model.O.density   = nO
    model.NOp.density = nNOp
    model.O2p.density = nO2p
    model.Op.density  = nOp
    
    model.Np.density  = model.N2.density*0
    model.N2p.density = model.N2.density*0
    model.N.density   = model.N2.density*0
    model.NO.density  = model.N2.density*0
    
    #model.e.density = np.sum([nNOp,nO2p,nOp,model.N2p.density], axis = 0)
    
    model.e.density = con["ne"][:, 0]
    
    model.check_chargeNeutrality()
    
    #assign production (unused?)
    for c in model.all_species:
        c.prod = e_prod *0
    
    #model.e.prod   = e_prod
    Op_prod  = e_prod * 0.56 * model.O.density  / (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)
    O2p_prod = e_prod * 1.00 * model.O2.density / (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)
    N2p_prod = e_prod * 0.92 * model.N2.density / (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)
    
    
    #create smooth function for production
    t = np.arange(0, te[-1], 0.01)
    
    stepf = np.array([stepped_prod_t(e_prod, i) for i in t])
    e_prod_smooth = PchipInterpolator(t, stepf)
    
    stepf = np.array([stepped_prod_t(Op_prod, i) for i in t])
    Op_prod_smooth = PchipInterpolator(t, stepf)
    
    stepf = np.array([stepped_prod_t(O2p_prod, i) for i in t])
    O2p_prod_smooth = PchipInterpolator(t, stepf)

    stepf = np.array([stepped_prod_t(N2p_prod, i) for i in t])
    N2p_prod_smooth = PchipInterpolator(t, stepf)
    
    #zero function for species that dont feel production 
    def zerof(t):
        try:
            some_object_iterator = iter(t)
            return np.zeros(len(t), len(z_model))
        except TypeError as te:
            return np.zeros(len(z_model))
    
    #save production function as array, in order of species
    prodMat = np.array([e_prod_smooth, 
                        zerof, 
                        Op_prod_smooth, 
                        zerof, 
                        O2p_prod_smooth, 
                        zerof, 
                        zerof, 
                        zerof, 
                        N2p_prod_smooth, 
                        zerof, 
                        zerof])

    #new way to do ionospheric chemistry
    #1. take all info from reaction (reaction rate, constituents)
    ode_mat = np.zeros((len(model.all_reactions), len(model.all_species)), dtype = 'object')
    rrate = np.array([np.array([r.r_rate_t(*t) for r in model.all_reactions]) for t in temp])
    
    for r in model.all_reactions:
        e = r.educts_ID
        p = r.products_ID
        ode_mat[r.r_ID, e[0]] = np.array([r.r_ID, -1, *e])
        ode_mat[r.r_ID, e[1]] = np.array([r.r_ID, -1, *e])
        ode_mat[r.r_ID, p[0]] = np.array([r.r_ID,  1, *e])
        ode_mat[r.r_ID, p[1]] = np.array([r.r_ID,  1, *e])
        if p[1] == p[0]: ode_mat[r.r_ID, p[0]] = np.array([r.r_ID,  2, *e])
    
    #2. produce raw DG from 1., excluding all terms that are 0 anyways.
    ode_raw = np.empty(len(model.all_species), dtype = 'object')
    for i in model.all_species:
    #    print(np.array([o for o in ode_mat[:, i.c_ID] if type(o)!= int]))
        ode_raw[i.c_ID] = np.array([o for o in ode_mat[:, i.c_ID] if type(o) != int])
        
    #3. produce DG with only relevant terms
    def fun(t, n, h):
        if t >= ts[0] and t<=te[-1]:
            k = len(ts[ts<=t])-1
        else:
            print(t)
            raise RuntimeError
            
        dndt = np.array([np.sum(((rrate[k, ode_raw[i.c_ID][:, 0], h].T  *   ode_raw[i.c_ID][:, 1] ).T                                    * n[ode_raw[i.c_ID][:, 2]]          * n[ode_raw[i.c_ID][:, 3]]                                   ), axis = 0)                         + prodMat[i.c_ID](t)[h]                        for i in model.all_species])
        return dndt
    
    def ic():
        res = np.empty(model.n_heights, dtype = 'object')
        
        for h in range(model.n_heights):
            n = np.array([c.density[h] for c in model.all_species])
            res[h] = solve_ivp(fun, (ts[0], te[-1]), n, method='BDF',vectorized=False, args = [h], t_eval = ts, max_step = 0.0444)
    #        res[h] = solve_ivp(fun, (ts[0], te[-1]), n, method='BDF',vectorized=False, args = [h], t_eval = np.arange(0, te[-1], 0.01), max_step = 0.0444)
            #for j, c in enumerate(model.all_species):
            #    c.density[h] = res[h].y[j, -1]           #replaced by below code:
        
        new_densities = np.array([alt.y[:, -1] for alt in res]).T
        #return new_densities
        return res
    
    res = ic()

    for h, i in enumerate(res):
        for c in model.all_species:
            plt.figure()
            plt.plot(i.t, i.y[c.c_ID, :], label = c.name)
            if c == model.e: plt.plot(ts, con["ne"][h, :], label = 'ElSpec ne')
            if c == model.N2: plt.plot(ts, con["iri"][h, 3], label = 'ElSpec N2')
            plt.legend(loc = 2)
            plt.yscale('log')
            ax2 = plt.gca().twinx()
            ax2.plot(ts, e_prod[:, 0], '.', color = 'orange', label = 'q_e')
            ax2.set_yscale('log')
            ax2.legend(loc = 1)
            #for t in ts: plt.axvline(t, alpha = 0.1)
        break
        
    elspec_iri_ICsorted = np.array([r.y for r in res])

    eff_rr = (rrate.T[:, 0, :]*elspec_iri_ICsorted[:, 10, :] +               rrate.T[:, 1, :]*elspec_iri_ICsorted[:, 4 , :] +               rrate.T[:, 2, :]*elspec_iri_ICsorted[:, 6 , :]   ) / elspec_iri_ICsorted[:, 0, :]

    mdict = {"elspec_iri_ICsorted": elspec_iri_ICsorted, "eff_rr": eff_rr}
    spio.savemat('IC.mat', mdict)


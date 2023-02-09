#!/usr/bin/env python
# coding: utf-8

# In[2]:


import scipy.io
import numpy as np
from os import listdir
from os.path import isfile, join


def loadEISCATdata(directory, identifiers):
    """
    loads a list of variables from a directorz containing EISCAT data files
    Parameters:
    directory: path to the EISCAT .mat files
        string
    identifiers: Variables as a list of strings
        ndarray
    """
    data = np.empty(len(identifiers), dtype = 'object')
    for i, var in enumerate(identifiers):
        data_files = [f for f in listdir(directory) if isfile(join(directory, f)) and f[-3:] == 'mat']
        n = len(data_files)
        mat_data = [scipy.io.loadmat(directory + file) for file in sorted(data_files)]
        data[i] = mat_extract_variable(mat_data, var)
        del mat_data
    return data


def mat_extract_variable(loaded_mat, var):
    """
    Extracts data for variable var from pre-loaded matfiles:
    loaded_mat = scipy.io.loadmat(path_to_file)
    data = loaded_mat('var')
    CAN HANDLE JAGGED ARRAYS (at least in 1 dimension)
    
    Parameters:
    loaded_mat: preloaded mat files as a list
        ndarry
    var: Variable to be extracted
        string
        
    Returns:
    data: all data saved in the loaded_mat['var']
    """
    n = len(loaded_mat)
    var_data = np.array([data[var] for data in loaded_mat])
    
    max_shape = np.amax([i.shape for i in var_data], axis = 0)
    array_size = [n, *max_shape]

    data = np.empty(array_size, dtype = 'object')
    
    for i in range(n):
        data[i, :len(var_data[i])] = var_data[i]
    
    if max_shape.shape == (2,):
        if max_shape[1] == 1:
            data = data.reshape(array_size[:-1])
    return data.astype('float')


# In[5]:





# In[ ]:





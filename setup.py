#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from datetime import datetime
from shutil import copy
import pickle


class print_config():
    """
    Controls print settings. This module is imported in all others (that can print),
    and share one instance of this class, defined in the main module.
    """
    def __init__(self, printing = 0):
        self.printing = printing
        self.print = self.printing
        
    def enable(self):
        self.printing = 1
        self.print = self.printing
    
    def disable(self):
        self.printing = 0
        self.print = self.printing
        

pcfg = print_config()


class setup:
    def __init__(self, msis_config, iri_config, chemistry_config, path_eiscat_data):
        #self._directory = os.path.abspath('')
        self._create_directory()
        self._copy_config(msis_config, iri_config, chemistry_config, path_eiscat_data)
     
    def _create_directory(self):
        """
        Creates a directory to store the configuration and results of this run.
        The directory is labeled with the time of starting the program.
        Folders created are Simulation (Config is created at the time of copying the config files)
        """
        self._today = datetime.now()
        #self._log_directory = self._directory + '/log/testing/' + self._today.strftime('%Y.%m.%d_%H_%M_%S')
        self._log_directory = 'log/testing/' + self._today.strftime('%Y.%m.%d_%H_%M_%S')
        if not os.path.isdir('log/testing/'):
            os.mkdir('log/testing/')
        os.mkdir(self._log_directory)
        os.mkdir(self._log_directory + '/Simulation')
        os.mkdir(self._log_directory + '/plots')

    def _copy_config(self, msis_config, iri_config, chemistry_config, path_eiscat_data):
        """
        Copies all configuration files into the log folder
        """
        os.mkdir(self._log_directory + '/Config')
        copy(msis_config, self._log_directory + '/Config')
        copy(iri_config, self._log_directory + '/Config')
        copy(chemistry_config, self._log_directory + '/Config')
        with open(self._log_directory + '/Config/path_to_eiscat_data.txt', 'w') as f:
            f.write(path_eiscat_data)
            
    def datadump(self, data):
        '''
        Dumps log data into the log file
        '''
        with open(self._log_directory + '/Simulation/log.p', 'ab') as f:
            pickle.dump(data, f)


# In[2]:


""" testing
msis_config = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/msis.txt'
iri_config  = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/iri.txt'
chemistry_config = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Meta-data/Reaction rates.txt'
path_eiscat_data = '/Users/ost051/Documents/PhD/Electron Precipitation/example/Data/'

setup_ = setup(msis_config, iri_config, chemistry_config, path_eiscat_data)

setup_.datadump(0, 1)
"""


# In[ ]:





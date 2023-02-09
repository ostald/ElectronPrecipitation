#!/usr/bin/env python
# coding: utf-8

# In[1]:


#toy model ODE


# In[2]:


from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'widget')


# In[3]:


def fun(t, y):
    return np.array([ -1*y[0]*y[1] + 2000*y[0]*y[2] ,
                     -2*y[1]*y[2],
                     -2000*y[2]*y[0]
                    ])
    #return np.array([-1, -2, -100]) * y


# In[4]:


res = solve_ivp(fun, (0, 10), np.array([1, 1, 1]), method = 'Radau')


# In[5]:


res


# In[6]:


res.t
res.y[:, -1]


# In[7]:


plt.figure()
plt.plot(res.t, res.y[0], label= '0')
plt.plot(res.t, res.y[1], label= '1')
plt.plot(res.t, res.y[2], label= '2')
plt.yscale('log')
plt.legend()


# In[8]:


1/np.exp(3)


# In[ ]:





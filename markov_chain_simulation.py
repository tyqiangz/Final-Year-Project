#!/usr/bin/env python
# coding: utf-8

# In[3]:


from FYP import *


# In[4]:


# Setting parameters ###################################

# code parameters
n0 = 2
r = 1019
wi = 39

#decryption parameters
t = 27

samp_method = 3

################################## Processing parameters ###################################

n = n0 * r
w = wi * n0


# In[5]:


print("(r,d,t) = (%d,%d,%d)" % (r,wi,t))
rhobar = rhoBar(n,w,t)
get_ipython().run_line_magic('time', 'prob = syndromeDist(n,w,t, rhobar)')


# In[6]:


t_pass = 10
t_fail = 30
print("############ DFR ################")
print("(t_pass, t_fail, r, d, t_init, samp_method) = (%d, %d, %d, %d, %d, %d)" % (t_pass, t_fail, r, wi, t, samp_method) )
get_ipython().run_line_magic('time', 'print("\\n\\nMC DFR =", DFR_model(t_pass, t_fail, n, w, t, prob, samp_method))')


# In[ ]:





# In[7]:


# Setting parameters ###################################

# code parameters
n0 = 2
r = 1171
wi = 39

#decryption parameters
t = 27

samp_method = 3

################################## Processing parameters ###################################

n = n0 * r
w = wi * n0


# In[8]:


print("(r,d,t) = (%d,%d,%d)" % (r,wi,t))
rhobar = rhoBar(n,w,t)
get_ipython().run_line_magic('time', 'prob = syndromeDist(n,w,t, rhobar)')


# In[9]:


t_pass = 10
t_fail = 40
print("############ DFR ################")
print("(t_pass, t_fail, r, d, t_init, samp_method) = (%d, %d, %d, %d, %d, %d)" % (t_pass, t_fail, r, wi, t, samp_method) )
get_ipython().run_line_magic('time', 'print("\\n\\nMC DFR =", DFR_model(t_pass, t_fail, n, w, t, prob, samp_method))')


# In[ ]:





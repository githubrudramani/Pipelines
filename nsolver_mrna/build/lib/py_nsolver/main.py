#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Using this guideline: https://nanostring.com/wp-content/uploads/Gene_Expression_Data_Analysis_Guidelines.pdf
# for geNorm (normalization )::?https://tep.cancer.illinois.edu/files/2020/08/MAN-10030-03_nCounter_Advanced_Analysis_2.0_User_Manual-1.pdf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.stats import linregress as  lr
import glob
import seaborn as sns
from scipy.stats import ttest_ind, ttest_rel
from statsmodels.stats.multitest import fdrcorrection

# import functions from utility module
import utility


# In[4]:


# Let's see what is in utility
help(utility)


# In[5]:


# first read the RCC file from read_rcc function from utility module
# let's see required arguments for the function
print(utility.read_rcc.__doc__)


# You can minimize the number of false positives by setting the threshold to the mean plus two
# standard deviations or to the maximum of the negative control counts; note, however, that
# although false positives will be rare, false negatives may be relatively abundant.
# 
# Conversely, you can set a more liberal threshold, such as the geometric mean of the negative
# controls. This will increase the number of false positives, butsimultaneously decrease the number
# of false negatives. 

# In[6]:


path = "../../RE__Murine_nanostring_data_from_sulforaphane_experiments/"

# It has two argumens 1. path:  path of RCC files and 2. bc : background correction method
raw_count, bc_count,  QC, PC = utility.read_rcc(path = path, bc = "mean")


# In[7]:


raw_count.head()


# In[8]:


bc_count.head()


# In[9]:


QC.head()


# In[10]:


# check the concentration for POS_E (0.5) and POS_F(0.125) values are correct for not
PC.loc[["POS_E", "POS_F"]].head()


# In[11]:


PC.loc[["POS_E", "POS_F"]].tail()


# In[12]:


# Lets inspect the sample using sample_qc function
print(utility.sample_qc.__doc__)


# In[13]:


QC_filtered, QC_unfiltered = utility.sample_qc(QC, fov_min = 95)


# In[14]:


QC_filtered.head()


# In[15]:


# Let's Normalize the count data
print(utility.normalize_data.__doc__)


# In[17]:


# If you want to remove the samples that does not pass QC uncomment below code
# sample = list(QC_filtered.index)
# sample.insert(0,"class")
# bc_count_filtered = bc_count[sample]
# bc_count_filtered
# Use bc_count_filtered instead of bc_count data

count_norm = utility.normalize_data(data = bc_count, method = "geNorm")


# In[18]:


# Let's perform two sample t-test
print(utility.DE_analysis.__doc__)


# In[20]:


# Let's load the sample metadata
sample = pd.read_csv("../sample.csv", index_col = 0)
sample.head()


# In[23]:


stat = utility.DE_analysis(count_norm, sample, sample_col = "treatment", group1 = "NS", group0 = "NV", ttest = "ind")
#stat.to_csv("stats.csv")
stat


# In[24]:


# Let's plot volcano
print(utility.plot_volcano.__doc__)


# In[25]:


utility.plot_volcano(stat)


# In[20]:


# Let's see up regulated and down regulated genes
stat = stat.sort_values("lfc", ascending = False)
stat_top = stat[stat.p_value < 0.05]
#stat_top.to_csv("top_genes.csv")
stat_top


# ### For Pathways analysis use my script in R studio at:
# ####  https://github.com/githubrudramani/Pipelines/blob/main/Bulk-RNAseq/PathwaysAnalysis/Enrichment.R

# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# ### Commonly used auxiliary functions

# In[1]:


import json

def isFloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

def loadJSONFile(path_to_file):
    with open(path_to_file) as f:
        try:
            loaded_json = json.load(f)
            return loaded_json
        except ValueError:
            return None

def writeJSONFile(data_dictionary, path_to_file):
    with open(path_to_file, 'w') as f:
            json.dump(data_dictionary, f)


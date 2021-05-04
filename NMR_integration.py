#!/usr/bin/env python
# coding: utf-8

# # Processing kinetic NMR integration data
# 
# 0. Imports

# In[1]:


import csv


# 1. The following is an example of the NMR integration data csv file. The second column of the csv file is the time measurement in seconds. The third column is the NMR integral value.

# ,,,,,,
# 
# #,X(I),Y(X),Y'(X),,,
# 
# Model,ARR_DATA(I),"Integral(10.013,9.493)",,,,
# 
# 1,0.000,-457.812,0.000,,,
# 
# 2,225.000,275.346,0.000,,,
# 
# 3,484.000,392.835,0.000,,,
# 
# 4,663.000,567.457,0.000,,,
# 
# 5,842.000,768.767,0.000,,,
# 
# 6,1021.000,1119.643,0.000,,,
# 
# 7,1201.000,1269.151,0.000,,,
# 
# 8,1380.000,1562.304,0.000,,,
# 
# 9,1559.000,1799.376,0.000,,,
# 
# 10,1738.000,2261.593,0.000,,,
# 
# 11,1917.000,2504.195,0.000,,,

# In[1]:


def readNMRIntegrationFromFile(path_to_file):
    '''Extracts the NMR integration time-series data from the csv file

    Parameters
    ----------
    path_to_file : string
        Path used to locate the NMR integration csv file.

    Returns
    -------
    array of tuples
       NMR integration time-series data
       The first element of each tuple is the time point in seconds. 
       The second element is the NMR integral value.
    '''
    nmr_integrals = []
    skip_number_lines = 8
    with open(path_to_file, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if (reader.line_num > skip_number_lines) and (reader.line_num % 2 == 1):
                nmr_integrals.append((float(row[1]), float(row[2])))
    return nmr_integrals


# In[ ]:


def readNMRSpectraFromFile(path_to_file):
    chemical_shifts = []
    nmr_spectra = []
    skip_number_lines = 1
    with open(path_to_file, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if (reader.line_num > skip_number_lines):
                chemical_shifts.append(float(row[0]))
                nmr_spectra.append([float(intensity) for intensity in row[1:-1]])
    return chemical_shifts, nmr_spectra


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# # Reading and writing structured data from experiments as JSON
# 
# 0. Imports

# In[1]:


import os
from collections import defaultdict
from .util import loadJSONFile, writeJSONFile


# ## Experiments: 
# #### I. Determining concentrations by absorption spectroscopy using principal component analysis on Fourier-transformed spectra.

# In[2]:


def standardSpectrumToJSON(component_ids, concentrations, spectrum):
    '''Structures the data obtained after measuring the spectrum of a standard mixture.

    Parameters
    ----------
    component_ids : array of strings
        Identifier associated to each component involved in the calibration experiment.
    concentrations : array of floats
        Concentration of each component in units of moles per liter. The concentration value of a
        component must be in the same array position as the component identifier in the component_ids array.
    spectrum : array of tuples
        Absorption spectrum obtained from the standard arrayed as follows.
        [(wavelength_1, absorbance_1), (wavelength_2, absorbance_2), ...]
        Wavelength and absorbance values are floats.
    Returns
    -------
    dictionary
       input data structured as a dictionary that can be readily written as a JSON file.
       Example:
       {
        "components": [{ id: component_1, concentration: conc_1 }, ...],
        "wavelengths": [wavelength_1, wavelength_2, ..., wavelength_i],
        "absorbances": [abs_1, abs_2, ..., abs_i]
       }
    '''
    assert (len(component_ids) == len(concentrations)),'''
    There should be a one-to-one correspondence between components and concentrations.
    '''
    wavelengths, absorbances = list(zip(*spectrum))
    json_dictionary = {
        "components": list(map(
            lambda component_id, concentration:
            {
                "id": component_id,
                "concentration": concentration
            },
            component_ids,
            concentrations
        )),
        "wavelengths": wavelengths,
        "absorbances": absorbances
    }
    return json_dictionary


# In[3]:


def readCalibrationData(path_to_directory):
    '''Reads the JSON files containing spectra obtained during the calibration experiments
       and merges the into single dictionary structure.

    Parameters
    ----------
    path_to_directory : string
        Path used to locate the directory containing the calibration .json files.
    Returns
    -------
    dictionary
       merged calibration data
       Example:
       {
        "component_ids": [component_1, component_2, ..., component_m],
                         // standard_1, standard_2, ..., standard_n  
        "concentrations": [[conc_1_1,   conc_1_2,   ..., conc_1_n], // component_1
                           [conc_2_1,   conc_2_2,   ..., conc_2_n], // component_2
                                                    ...,
                           [conc_m_1,   conc_m_2,   ..., conc_m_n], // component_m
        "wavelengths": [wavelength_1, wavelength_2, ..., wavelength_i],
                      // standard_1,  standard_2, ...,  standard_n  
        "absorbances": [[abs_1_1,     abs_1_2,    ...,  abs_1_n], // wavelength_1
                        [abs_2_1,     abs_2_2,    ...,  abs_2_n], // wavelength_2
                                                  ...,
                        [abs_m_1,     abs_m_2,    ...,  abs_m_n], // wavelength_i  
        }
    '''
    # List the standard spectrum JSON files in the directory, ignoring hidden files
    calibration_spectrum_files = [f for f in os.listdir(path_to_directory) if not f.startswith('.')]
    merged_spectra = defaultdict(list)
    merged_concentrations = defaultdict(list)
    
    # Iterate over each standard spectrum file in the directory.
    for file_name in calibration_spectrum_files:
        path_to_file = os.path.join(path_to_directory, file_name)
        standard_spectrum_data = loadJSONFile(path_to_file)
        if standard_spectrum_data is not None:
            spectrum = list(zip(standard_spectrum_data["wavelengths"], standard_spectrum_data["absorbances"]))
            for wavelength, absorbance in spectrum:
                merged_spectra[wavelength].append(absorbance)
            components = list(map(lambda component: tuple(component.values()), standard_spectrum_data["components"]))
            for component_id, concentration in components:
                merged_concentrations[component_id].append(concentration)
        
    # Check that all components have a defined concentration value in each standard spectrum file
    number_of_standards = len(calibration_spectrum_files)
    standards_with_component = iter(merged_concentrations.values())
    if not all(len(standards) == number_of_standards for standards in standards_with_component):
         raise ValueError('All components must have a defined concentration value in each standard spectrum file')

    component_ids, concentrations =  zip(*merged_concentrations.items())
    
    # Keep only the absorbances measured at the intersection of the set of wavelengths across all standards
    wavelengths, absorbances = zip(*filter(
        lambda absorbances_at_wavelength:
        len(absorbances_at_wavelength) != number_of_standards,
        sorted(merged_spectra.items())
    ))
    
    return {
        "component_ids": list(component_ids),
        "concentrations": list(concentrations),
        "wavelengths": list(wavelengths),
        "absorbances": list(absorbances)
    }


# In[4]:


def kineticSpectrumToJSON(experiment_id, time_point, dilution_factor, spectrum):
    '''Structures the data obtained after measuring the spectrum of an unknown mixture in a kinetic experiment.

    Parameters
    ----------
    experiment_id : string
        Identifier associated to the kinetic experiment.
    time_point : float
        Time measurement in seconds taken during the experiment.
    dilution_factor: float
        Dilution factor = (total volume in cuvette)/(sample aliquot volume)
    spectrum : array of tuples
        Absorption spectrum obtained from the unknown mixture arrayed as follows.
        [(wavelength_1, absorbance_1), (wavelength_2, absorbance_2), ...]
        Wavelength and absorbance values are floats.
    Returns
    -------
    dictionary
       input data structured as a dictionary that can be readily written as a JSON file.
       Example:
       {
        "experiment_id": X_1234,
        "time_point": 60.0,
        "dilution_factor": 3.0,
        "wavelengths": [wavelength_1, wavelength_2, ..., wavelength_i],
        "absorbances": [abs_1, abs_2, ..., abs_i]
       }
    '''
    wavelengths, absorbances = list(zip(*spectrum))
    json_dictionary = {
        "experiment_id": experiment_id,
        "time_point": time_point,
        "dilution_factor": dilution_factor,
        "wavelengths": wavelengths,
        "absorbances": absorbances
    }
    return json_dictionary


# In[1]:


def readKineticData(path_to_directory):
    '''Reads the JSON files containing spectra obtained during the kinetic experiments
       and merges the into single dictionary structure.

    Parameters
    ----------
    path_to_directory : string
        Path used to locate the directory containing the kinetic .json files.
    Returns
    -------
    dictionary
       merged kinetic data
       Example:
       {
        "experiment_id": X_1234,
        "time_points": [t_1, t_2, ..., t_n],
        "dilution_factors": [d_1, d_2, ..., d_n],
        "wavelengths": [wavelength_1, wavelength_2, ..., wavelength_i],
                      // at_t_1,    at_t_2,    ..., at_t_n  
        "absorbances": [[abs_1_1,   abs_1_2,   ..., abs_1_n], // wavelength_1
                        [abs_2_1,   abs_2_2,   ..., abs_2_n], // wavelength_2
                                               ...,
                        [abs_m_1,   abs_m_2,   ..., abs_m_n], // wavelength_i  
        }
    '''
    # List the spectrum files in the directory associated to the experiment, ignoring hidden files
    kinetic_spectrum_files = [f for f in os.listdir(path_to_directory) if not f.startswith('.')]
    
    merged_spectra = defaultdict(list)
    time_points = []
    experiment_ids = []
    dilution_factors = []
    
    # Iterate over each spectrum file in the directory.
    for file_name in kinetic_spectrum_files:
        path_to_file = os.path.join(path_to_directory, file_name)
        kinetic_spectrum_data = loadJSONFile(path_to_file)
        if kinetic_spectrum_data is not None:
            spectrum = list(zip(kinetic_spectrum_data["wavelengths"], kinetic_spectrum_data["absorbances"]))
            for wavelength, absorbance in spectrum:
                merged_spectra[wavelength].append(absorbance)
            time_points.append(kinetic_spectrum_data["time_point"])
            experiment_ids.append(kinetic_spectrum_data["experiment_id"])
            dilution_factors.append(kinetic_spectrum_data["dilution_factor"])
    
    indexes = [i for (v, i) in sorted((v, i) for (i, v) in enumerate(time_points))]
        
    # Check that all experiments have the same id in each standard spectrum file
    number_of_experiments = len(experiment_ids)
    id_iter = iter(experiment_ids)
    experiment_id = next(id_iter)
    if not all(identifier == experiment_id for identifier in id_iter):
         raise ValueError('All experiments must have the same id.')
            
    # Sort time points in ascending order and the corresponding kinetic data accordingly
    time_points = [time_points[i] for i in indexes]
    dilution_factors = [dilution_factors[i] for i in indexes]
    
    # Keep only the absorbances measured at the intersection of the set of wavelengths across all standards
    wavelengths, absorbances = zip(*map(
        lambda absorbances_at_wavelength:
        (absorbances_at_wavelength[0], [absorbances_at_wavelength[1][i] for i in indexes]),
        filter(
            lambda absorbances_at_wavelength:
            len(absorbances_at_wavelength) != number_of_experiments,
            sorted(merged_spectra.items())
        )))
    
    
    return {
        "experiment_id": experiment_id,
        "time_points": time_points,
        "wavelengths": list(wavelengths),
        "absorbances": list(absorbances),
        "dilution_factors": dilution_factors
    }


# #### II. Determining activation energy from NMR integration time-series data at different temperatures.

# In[2]:


def NMRIntegrationToJSON(experiment_id, temperature, nmr_integrals):
    '''Structures the NMR integration data obtained from a time-series experiment carried out at a certain temperature.

    Parameters
    ----------
    experiment_id : string
        Identifier associated to the kinetic experiment.
    temperature : float
        Temperature measurement in kelvins.
    nmr_integrals : array of tuples
        NMR integral values arrayed as follows.
        [(t_1, integral_1), (t_2, integral_2), ...]
        Time and integral values are floats.
    Returns
    -------
    dictionary
       input data structured as a dictionary that can be readily written as a JSON file.
       Example:
       {
        "experiment_id": X_1234,
        "temperature": 250.0,
        "time_points": [t_1, t_2, ..., t_i],
        "integrals": [integral_1, integral_2, ..., integral_i]
       }
    '''
    time_points, integrals = list(zip(*nmr_integrals))
    json_dictionary = {
        "experiment_id": experiment_id,
        "temperature": temperature,
        "time_points": time_points,
        "integrals": integrals
    }
    return json_dictionary


# In[1]:


def readNMRIntegrationData(path_to_file):
    '''Reads the JSON file containing NMR integration data obtained during a time-series experiment
    carried out at a certain temperature.

    Parameters
    ----------
    path_to_directory : string
        Path used to locate the directory containing the nmr integration .json files.
    Returns
    -------
    dictionary
       NMR integration data
       Example:
       {
        "experiment_id": X_1234,
        "temperature": 250.0,
        "time_points": [t_1, t_2, ..., t_i],
        "integrals": [integral_1, integral_2, ..., integral_i]
       }
    '''
    return loadJSONFile(path_to_file)


# In[2]:


def standardNMRSpectrumToJSON(component_ids, concentrations, spectrum):
    '''Structures the data obtained after measuring the NMR spectrum of a standard mixture.

    Parameters
    ----------
    component_ids : array of strings
        Identifier associated to each component involved in the calibration experiment.
    concentrations : array of floats
        Concentration of each component in units of moles per liter. The concentration value of a
        component must be in the same array position as the component identifier in the component_ids array.
    spectrum : array of tuples
        NMR spectrum obtained from the standard arrayed as follows.
        [(shift_1, intensity_1), (shift_2, intensity_2), ...]
        Chemical shift and intensity values are floats.
    Returns
    -------
    dictionary
       input data structured as a dictionary that can be readily written as a JSON file.
       Example:
       {
        "components": [{ id: component_1, concentration: conc_1 }, ...],
        "chemical_shifts": [shift_1, shift_2, ..., shift_i],
        "intensities": [intensity_1, intensity_2, ..., intensity_i]
       }
    '''
    assert (len(component_ids) == len(concentrations)),'''
    There should be a one-to-one correspondence between components and concentrations.
    '''
    chemical_shifts, intensities = list(zip(*spectrum))
    json_dictionary = {
        "components": list(map(
            lambda component_id, concentration:
            {
                "id": component_id,
                "concentration": concentration
            },
            component_ids,
            concentrations
        )),
        "chemical_shifts": chemical_shifts,
        "intensities": intensities
    }
    return json_dictionary


# In[3]:


def kineticNMRSpectrumToJSON(experiment_id, time_point, dilution_factor, spectrum):
    '''Structures the data obtained after measuring the spectrum of an unknown mixture in a kinetic experiment.

    Parameters
    ----------
    experiment_id : string
        Identifier associated to the kinetic experiment.
    time_point : float
        Time measurement in seconds taken during the experiment.
    dilution_factor: float
    spectrum : array of tuples
        NMR spectrum obtained from the unknown mixture arrayed as follows.
        [(shift_1, intensity_1), (shift_2, intensity_2), ...]
        Chemical shift and intensity values are floats.
    Returns
    -------
    dictionary
       input data structured as a dictionary that can be readily written as a JSON file.
       Example:
       {
        "experiment_id": X_1234,
        "time_point": 60.0,
        "dilution_factor": 3.0,
        "chemical_shifts": [shift_1, shift_2, ..., shift_i],
        "intensities": [intensity_1, intensity_2, ..., intensity_i]
       }
    '''
    chemical_shifts, intensities = list(zip(*spectrum))
    json_dictionary = {
        "experiment_id": experiment_id,
        "time_point": time_point,
        "dilution_factor": dilution_factor,
        "chemical_shifts": chemical_shifts,
        "intensities": intensities
    }
    return json_dictionary


# In[4]:


def readNMRCalibrationData(path_to_directory):
    '''Reads the JSON files containing spectra obtained during the calibration experiments
       and merges the into single dictionary structure.

    Parameters
    ----------
    path_to_directory : string
        Path used to locate the directory containing the calibration .json files.
    Returns
    -------
    dictionary
       merged calibration data
       Example:
       {
        "component_ids": [component_1, component_2, ..., component_m],
                         // standard_1, standard_2, ..., standard_n  
        "concentrations": [[conc_1_1,   conc_1_2,   ..., conc_1_n], // component_1
                           [conc_2_1,   conc_2_2,   ..., conc_2_n], // component_2
                                                    ...,
                           [conc_m_1,   conc_m_2,   ..., conc_m_n], // component_m
        "chemical_shifts": [shift_1, shift_2, ..., shift_i],
                      // standard_1,  standard_2, ...,  standard_n  
        "intensities": [[i_1_1,     i_1_2,    ...,  i_1_n], // shift_1
                        [i_2_1,     i_2_2,    ...,  i_2_n], // shift_2
                                                  ...,
                        [i_m_1,     i_m_2,    ...,  i_m_n], // shift_i  
        }
    '''
    # List the standard spectrum JSON files in the directory, ignoring hidden files
    calibration_spectrum_files = [f for f in os.listdir(path_to_directory) if not f.startswith('.')]
    merged_spectra = defaultdict(list)
    merged_concentrations = defaultdict(list)
    
    # Iterate over each standard spectrum file in the directory.
    for file_name in calibration_spectrum_files:
        path_to_file = os.path.join(path_to_directory, file_name)
        standard_spectrum_data = loadJSONFile(path_to_file)
        if standard_spectrum_data is not None:
            spectrum = list(zip(standard_spectrum_data["chemical_shifts"], standard_spectrum_data["intensities"]))
            for chemical_shift, intensity in spectrum:
                merged_spectra[chemical_shift].append(intensity)
            components = list(map(lambda component: tuple(component.values()), standard_spectrum_data["components"]))
            for component_id, concentration in components:
                merged_concentrations[component_id].append(concentration)
        
    # Check that all components have a defined concentration value in each standard spectrum file
    number_of_standards = len(calibration_spectrum_files)
    standards_with_component = iter(merged_concentrations.values())
    if not all(len(standards) == number_of_standards for standards in standards_with_component):
         raise ValueError('All components must have a defined concentration value in each standard spectrum file')

    component_ids, concentrations =  zip(*merged_concentrations.items())
    
    # Keep only the absorbances measured at the intersection of the set of wavelengths across all standards
    chemical_shifts, intensities = zip(*filter(
        lambda intensities_at_shift:
        len(intensities_at_shift) != number_of_standards,
        sorted(merged_spectra.items())
    ))
    
    return {
        "component_ids": list(component_ids),
        "concentrations": list(concentrations),
        "chemical_shifts": list(chemical_shifts),
        "intensities": list(intensities)
    }


# In[ ]:


def readNMRKineticData(path_to_directory):
    '''Reads the JSON files containing spectra obtained during the kinetic experiments
       and merges the into single dictionary structure.

    Parameters
    ----------
    path_to_directory : string
        Path used to locate the directory containing the kinetic .json files.
    Returns
    -------
    dictionary
       merged kinetic data
       Example:
       {
        "experiment_id": X_1234,
        "time_points": [t_1, t_2, ..., t_n],
        "dilution_factors": [d_1, d_2, ..., d_n],
        "chemical_shifts": [shift_1, shift_2, ..., shift_i],
                      // standard_1,  standard_2, ...,  standard_n  
        "intensities": [[i_1_1,     i_1_2,    ...,  i_1_n], // shift_1
                        [i_2_1,     i_2_2,    ...,  i_2_n], // shift_2
                                                  ...,
                        [i_m_1,     i_m_2,    ...,  i_m_n], // shift_i
        }
    '''
    # List the spectrum files in the directory associated to the experiment, ignoring hidden files
    kinetic_spectrum_files = [f for f in os.listdir(path_to_directory) if not f.startswith('.')]
    
    merged_spectra = defaultdict(list)
    time_points = []
    experiment_ids = []
    dilution_factors = []
    
    # Iterate over each spectrum file in the directory.
    for file_name in kinetic_spectrum_files:
        path_to_file = os.path.join(path_to_directory, file_name)
        kinetic_spectrum_data = loadJSONFile(path_to_file)
        if kinetic_spectrum_data is not None:
            spectrum = list(zip(kinetic_spectrum_data["chemical_shifts"], kinetic_spectrum_data["intensities"]))
            for chemical_shift, intensity in spectrum:
                merged_spectra[chemical_shift].append(intensity)
            time_points.append(kinetic_spectrum_data["time_point"])
            experiment_ids.append(kinetic_spectrum_data["experiment_id"])
            dilution_factors.append(kinetic_spectrum_data["dilution_factor"])
    
    indexes = [i for (v, i) in sorted((v, i) for (i, v) in enumerate(time_points))]
        
    # Check that all experiments have the same id in each standard spectrum file
    number_of_experiments = len(experiment_ids)
    id_iter = iter(experiment_ids)
    experiment_id = next(id_iter)
    print(experiment_id)
    if not all(identifier == experiment_id for identifier in id_iter):
         raise ValueError('All experiments must have the same id.')
            
    # Sort time points in ascending order and the corresponding kinetic data accordingly
    time_points = [time_points[i] for i in indexes]
    dilution_factors = [dilution_factors[i] for i in indexes]
    
    # Keep only the absorbances measured at the intersection of the set of wavelengths across all standards
    chemical_shifts, intensities = zip(*map(
        lambda intensities_at_shift:
        (intensities_at_shift[0], [intensities_at_shift[1][i] for i in indexes]),
        filter(
            lambda intensities_at_shift:
            len(intensities_at_shift) != number_of_experiments,
            sorted(merged_spectra.items())
        )))
    
    
    return {
        "experiment_id": experiment_id,
        "time_points": time_points,
        "chemical_shifts": list(chemical_shifts),
        "intensities": list(intensities),
        "dilution_factors": dilution_factors
    }


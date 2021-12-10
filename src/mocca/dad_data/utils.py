#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 13:31:37 2021

@author: haascp
"""

import numpy as np

def get_compound_names(coumpound_data):
    # iterate through all values to find compounds with nonzero concs
    # input is a dictionary mapping compounds to concentrations
    return [compound for compound, conc in
            coumpound_data.compound_input.items() if conc != 0]

def get_compound_concentration(coumpound_data, compound_id):
    return coumpound_data.compound_input[compound_id]

def sum_absorbance_by_time(data):
    """
    Sums the absorbances for each time point over all wavelengths

    Parameters
    ----------
    data : numpy.ndarray
        Actual experimental data with shape [# of wavelengths] x [timepoints].
        Generated from dataframe with absorbance_to_array function

    Returns
    -------
    numpy.ndarray
        A 1D array containing the sum of wavelengths at each time point
    """
    return data.sum(axis=0)

def trim_data(data, time, length):
    """Trims the 2D DADData in the time dimension to the length provided"""
    if length < data.shape[1]:
        return data[:, :length], time[:length]
    else:
        return data, time

def absorbance_to_array(df):
    """
    Generates a 2D absorbance array of the absorbance values.

    Parameters
    ----------
    df : pandas.DataFrame
        Columns:
            time: Chromatogram time
            wavelength: Detection wavelength
            absorbance: absorbance value.

    Returns
    -------
    absorbance_array : numpy.ndarray
        Absorbance values with number of wavelengths in the first 
        and number of recorded times in the second dimension.

    """
    absorbance_array = df.absorbance.to_numpy().\
        reshape(df.wavelength.nunique(), df.time.nunique())
    return absorbance_array

def apply_filter(dataframe):
    """
    Filters absorbance data of tidy 3D DAD dataframes to remove noise
    and background systematic error.
    Parameters
    ----------
    dataframe : pandas.DataFrame
        Columns:
            time: Chromatogram time
            wavelength: Detection wavelength
            absorbance: absorbance value.

    Returns
    -------
    df : pandas.DataFrame
        Columns:
            time: Chromatogram time
            wavelength: Detection wavelength
            absorbance: absorbance value.
    """

    df = dataframe.copy()
    df['absorbance'] = df.groupby('time')['absorbance']\
        .rolling(window=5, center=True).mean().reset_index(0,drop=True)
    df = df.dropna().reset_index(0,drop=True)
    n_times = len(df.time.unique())
    wls = df.wavelength.unique()
    df['absorbance'] = df.absorbance - df[df['wavelength'] == wls.max()]\
        .absorbance.iloc[np.tile(np.arange(n_times), len(wls))]\
        .reset_index(0,drop=True)
    return df

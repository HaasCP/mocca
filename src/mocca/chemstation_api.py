#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 15:28:24 2021

@author: haascp
"""

# Files and folders navigation
import os

# Data manipulation
import numpy as np
import pandas as pd
import math
from micdrop.utils import read_csv

def read_csv_agilent(path):
    """
    Reads the UTF-16 encoded 3D data exported by the ChemStation macro.  
    Parameters
    ----------
    path : str
        The directory, in which the experimental data are stored.

    Returns
    -------
    df : pandas.DataFrame
        First column is time, the following columns obtain the absorbance 
        values at the given detection wavelength in the column name.
    """

    return read_csv(os.path.join(path, 'DAD1.CSV'), 'utf-16')

def tidy_csv_agilent(dataframe, wl_high_pass=None):
    """
    Tidies the raw data obtained from reading the CSV

    Parameters
    ----------
    dataframe : pandas.DataFrame
        First column is time, the following columns obtain the absorbance 
        values at the given detection wavelength in the column name.

    Raises
    ------
    ValueError
        If acquisition rate of the DAD was not constant, this error is raised.

    Returns
    -------
    df : pandas.DataFrame
        Columns:
            time: Chromatogram time
            wavelength: Detection wavelength
            absorbance: Absorbance value
    """

    df = dataframe.copy()
    df.rename(columns={df.columns[0]: 'time'}, inplace=True) # names time column
    
    acq_time = df.time.max() / len(df)
    
    time_series = pd.Series(range(1, (len(df) + 1))).astype(float) * acq_time # generates new time column
    df['time'] = time_series
    df = pd.melt(df, id_vars='time', value_vars=df.columns[1:], 
                 var_name='wavelength', value_name='absorbance')
    df['wavelength'] = df['wavelength'].astype(float)
    if wl_high_pass:
        df = df[df.wavelength >= wl_high_pass]
    acq_rate = 1 / (acq_time * 60)
    df.attrs = {'acq_time': acq_time,
                'acq_rate': round(acq_rate, 2 - int(math.floor(math.log10(abs(acq_rate)))) - 1)} # round to two significant digits calculates acquisition rate in Hz and safes it as dataframe metadata
    return df

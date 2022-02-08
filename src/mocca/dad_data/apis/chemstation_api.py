#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 15:28:24 2021

@author: haascp
"""

# Files and folders navigation
import os

# Data manipulation
import pandas as pd


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
    with open(os.path.join(path, 'DAD1.CSV'), 'r', encoding='utf-16') as f:
        df = pd.read_csv(f)
    return df

def tidy_df_agilent(dataframe, wl_high_pass=None, wl_low_pass=None):
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
    return df

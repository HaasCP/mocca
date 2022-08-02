# flake8: noqa
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 15:17:53 2021

@author: haascp
"""
import pandas as pd

from mocca.dad_data.utils import df_to_array, apply_filter


def read_txt_shimadzu(path):
    """
    Reads the 3D data exported by the LabSolutions software. 
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
    with open(path) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

    #cut text file to the relevant data 
    pda_line = 0
    while "[PDA 3D]" not in lines[pda_line]:
        pda_line += 1

    start_line = pda_line
    while len(lines[start_line]) < 100:
        start_line += 1
    data = lines[start_line:]

    #split data using comma as separator
    data = [line.split('\t') for line in data]
    
    #create dataframe and tidy data
    df = pd.DataFrame(data).dropna()
    #first line as column names
    df.columns = df.iloc[0]
    df = df.drop(df.index[0]).reset_index(drop=True)
    #first column time
    df.rename(columns={df.columns[0]: 'time'}, inplace=True) # names time column
    df = df.astype('float')
    #set time vector
    acq_time = df.time.max() / len(df)
    time_series = pd.Series(range(1, (len(df) + 1))).astype(float) * acq_time # generates new time column
    df['time'] = time_series
    #tidy data
    df = pd.melt(df, id_vars='time', value_vars=df.columns[1:], 
                 var_name='wavelength', value_name='absorbance')
    df['wavelength'] = df['wavelength'].astype(float)
    df['wavelength'] = df['wavelength'] / 100
    df['absorbance'] = df['absorbance'] / 1000
    return df


def read_labsolutions(path, wl_high_pass=None, wl_low_pass=None):
    """
    Labsolutions read and processing function.
    """
    df = read_txt_shimadzu(path)
    df = apply_filter(df, wl_high_pass, wl_low_pass)
    data, time, wavelength = df_to_array(df)
    return data, time, wavelength
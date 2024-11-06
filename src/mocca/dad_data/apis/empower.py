#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 14:53:54 2022

@author: haascp

For more information on how to export raw data out of Empower see
https://support.waters.com/KB_Inf/Empower_Breeze/WKB77571_How_to_export_3D_raw\
    _data_from_Empower_to_a_Microsoft_Excel_spreadsheet
"""

import numpy as np
import pandas as pd

from mocca.dad_data.utils import df_to_array, apply_filter


def read_arw_empower(path: str, wl_high_pass:float=None, wl_low_pass:float=None) -> pd.DataFrame:
    '''
    Read Empower ARW file and return pandas DataFrame.
    '''
    wavelength_vec = None
    # find line number where data starts and read wavelength vector
    with open(path) as f:
        for i, line in enumerate(f):
            if line.startswith('Wavelength'):
                wavelength_vec = line.split("\t")[1:]
            elif line.startswith('Time'):
                n_skip = i+1
                break
    # read data
    df = pd.read_csv(path, sep='\t', names=['time']+wavelength_vec, skiprows=n_skip, dtype=float)
    # melt data into 3 columns (time, wavelength, absorbance)
    df = pd.melt(df, id_vars='time', value_vars=df.columns,
                 var_name='wavelength', value_name='absorbance')
    return df


def read_empower(path, wl_high_pass=None, wl_low_pass=None):
    """
    Labsolutions read and processing function.
    """
    df = read_arw_empower(path)
    df = apply_filter(df, wl_high_pass, wl_low_pass)
    data, time, wavelength = df_to_array(df)
    return data, time, wavelength

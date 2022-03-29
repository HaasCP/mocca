#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 14:53:54 2022

@author: haascp
"""

"""
https://support.waters.com/KB_Inf/Empower_Breeze/WKB77571_How_to_export_3D_raw_data_from_Empower_to_a_Microsoft_Excel_spreadsheet

"""
import numpy as np
import pandas as pd

from mocca.dad_data.utils import df_to_array, apply_filter


def read_arw_empower(path, wl_high_pass=None, wl_low_pass=None):
    with open(path) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
    
    time_idx = [n for n, l in enumerate(lines) if l.startswith('Time')][0]
    wl_idx = [n for n, l in enumerate(lines) if l.startswith('Wavelength')][0]
    
    absorbance_list = []
    time_vec = []
    for line in lines[time_idx + 2:]:
        line_list = line.split("\t")
        absorbance_list.append(line_list[1:])
        time_vec.append(line_list[0])
        
    absorbance = np.array(absorbance_list).astype(float)
    time = [float(i) for i in time_vec]
    acq_time = max(time) / len(time)
    time_series = pd.Series(range(1, (len(time) + 1))).astype(float) * acq_time # generates new time column
    
    wavelength_vec = lines[wl_idx].split("\t")[1:]
    wavelength = [float(i) for i in wavelength_vec]
    
    df = pd.DataFrame(absorbance, columns=wavelength)
    df = df.assign(time=pd.Series(time_series).values)
    
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

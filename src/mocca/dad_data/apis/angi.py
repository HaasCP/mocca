#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/13/22

@author: haascp
"""
import pandas as pd

from mocca.dad_data.utils import df_to_array, apply_filter


def read_csv_angi(path):
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
    with open(path, 'r') as f:
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
    # name time column
    df.rename(columns={df.columns[0]: 'time'}, inplace=True)

    acq_time = df.time.max() / len(df)

    # generate new time column
    time_series = pd.Series(range(1, (len(df) + 1))).astype(float) * acq_time
    df['time'] = time_series
    df = pd.melt(df, id_vars='time', value_vars=df.columns[1:],
                 var_name='wavelength', value_name='absorbance')
    df['wavelength'] = df['wavelength'].astype(float)
    return df


def read_angi(path, wl_high_pass=None, wl_low_pass=None):
    """
    Chemstation read and processing function.
    """
    df = read_csv_angi(path)
    df = tidy_df_agilent(df)
    df = apply_filter(df, wl_high_pass, wl_low_pass)
    data, time, wavelength = df_to_array(df)
    return data, time, wavelength

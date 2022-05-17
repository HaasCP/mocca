# -*- coding: utf-8 -*-
"""
Created on Thu May 12 09:10:57 2022

@author: CPH
"""

import h5py
import numpy as np
import pandas as pd

from mocca.dad_data.utils import df_to_array, apply_filter


def read_adf_datacube(path):
    """
    Reads the raw data stored in the data cube layer, which are the HPLC-DAD
    absorbance values and the time scale.
    """
    f = h5py.File(path, 'r')
    data_cubes = f['data-cubes']
    uvvis_id = [idx for idx in data_cubes.keys() if "DAD1UV" in idx][0]
    uvvis_data = data_cubes[uvvis_id]
    measures = uvvis_data['measures']
    data_idx = list(measures.keys())[0]
    absorbance_idx = measures[data_idx]
    absorbance = absorbance_idx[()]
    absorbance = np.swapaxes(absorbance, 0, 1)

    scales = uvvis_data['scales']
    data_idx = list(scales.keys())[0]
    time_idx = scales[data_idx]
    time = time_idx[()]

    return np.array(absorbance), list(time)


def read_adf_description(wavelength_vals=None):
    """
    Queries the adf data description layer to extract the wavlength vector.
    For this query, the h5ld package is required which cannot be used on Windows
    machines currently. As a workaround, the user can give start and stop values
    as set on the DAD manually.
    """
    if wavelength_vals is not None:
        start, stop, length = wavelength_vals
        wavelength = list(np.linspace(start, stop, length))
    else:
        try:
            import h5ld  # noqa: F401
            # to be developed with Laura di Rocco
        except AttributeError:
            print("If the h5ld package cannot be installed on your machine, "
                  "you have to give the wavelength values manually for adf data.")
    return wavelength


def preprocess_df(df):
    """
    Preprocesses the df time column to be in line with the Chemstation API.
    """
    acq_time = df.time.max() / len(df)

    # generate new time column
    time_series = pd.Series(range(1, (len(df) + 1))).astype(float) * acq_time / 60
    df['time'] = time_series
    return df


def read_adf(path, wl_high_pass=None, wl_low_pass=None, wl_start=200, wl_stop=550):
    """
    Reads adf files as exported by the Agilent ADF Adapter.
    """
    absorbance, time = read_adf_datacube(path)
    wavelength_vals = (wl_start, wl_stop, absorbance.shape[0])
    wavelength = read_adf_description(wavelength_vals=wavelength_vals)

    df = pd.DataFrame(np.swapaxes(absorbance, 0, 1), columns=wavelength)
    df.insert(0, "time", time)
    df = preprocess_df(df)

    df = pd.melt(df, id_vars='time', value_vars=df.columns[1:],
                 var_name='wavelength', value_name='absorbance')
    df['wavelength'] = df['wavelength'].astype(float)
    df = apply_filter(df, wl_high_pass, wl_low_pass)
    data, time, wavelength = df_to_array(df)
    return data, time, wavelength

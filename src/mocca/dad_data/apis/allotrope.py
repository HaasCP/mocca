# -*- coding: utf-8 -*-
"""
Created on Thu May 12 09:10:57 2022

@author: CPH
"""

import numpy as np
import pandas as pd

from mocca.dad_data.utils import df_to_array, apply_filter


def get_uvvis_dataset_name(path):
    """
    Queries the data description layer of the adf file to find the name of the
    dataset which is of the type 'three-dimensional ultraviolet spectrum' as
    defined by the AFO.
    """
    from h5ld import AllotropeDF
    import h5py
    import rdflib
    
    with h5py.File(path, mode="r") as f:
        g = AllotropeDF(f).get_ld()
    
    datasets_query = '''SELECT ?s ?p ?o
                        WHERE { ?s ?p ?o .
                               FILTER regex(str(?o), "DataSet") .
                               }'''
    qres = g.query(datasets_query)

    datasets = [x[0] for x in qres]
    for d in datasets:
        subj = list(g.triples((None, None, d)))[1][0]
        dataset = list(
            g.triples(
            (subj, None, rdflib.term.URIRef('http://purl.allotrope.org/ontologies/result#AFR_0001527'))
            )
        )
        if dataset:
            dataset_name = d
    return dataset_name


def read_adf_datacube(path):
    """
    Reads the raw data stored in the data cube layer, which are the HPLC-DAD
    absorbance values and the time scale.
    """
    import h5py

    with h5py.File(path, mode="r") as f:
        data_cubes = f['data-cubes']
        dataset_name = get_uvvis_dataset_name(path)
        dataset_key = str(dataset_name).split(':')[1][2:].replace("/", '-')
        uvvis_data = data_cubes[dataset_key]
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
    For this query, the h5ld package is required which can be installed by 
    editable pip install from https://github.com/laura-dirocco/h5ld. In case there
    are problems with installation, the user can give start and stop values
    as set on the DAD manually.
    """
    if wavelength_vals is not None:
        start, stop, length = wavelength_vals
        wavelength = list(np.linspace(start, stop, length))
    else:
        try:
            pass
            # noqa: F401
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


def read_adf(path, wl_high_pass=None, wl_low_pass=None, wl_start=190, wl_stop=400):
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

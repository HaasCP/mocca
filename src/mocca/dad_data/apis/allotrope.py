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
            g.triples((subj, None,
                       rdflib.term.URIRef('http://purl.allotrope.org/ontologies/result#AFR_0001527')))  # noqa: E501
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


def get_function_paramenters(path):
    """
    Reads the parameters of the linear function which describes the wavelength
    vector out of the data description layer.
    """
    from h5ld import AllotropeDF
    import h5py
    import rdflib

    with h5py.File(path, mode="r") as f:
        g = AllotropeDF(f).get_ld()
    dataset = get_uvvis_dataset_name(path)

    # Step 1: Extract scaleMapping from a Dataset
    scale_mapping_id = list(g.objects(list(g.triples((None, None, dataset)))[0][0],
                                      rdflib.term.URIRef('http://purl.allotrope.org/ontologies/datacube-hdf-map#scaleMapping')))  # noqa: E501

    # In case of error, i.e., there is no FunctionScaleMapping, it returns 0, 0
    param1 = 0.0
    param2 = 0.0
    # Code based on the assumption that there is only one FunctionScaleMapping
    for mapping in scale_mapping_id:
        # Step 2: Looking for a Function Scale Mapping
        mapping_function = list(g.triples((mapping,
                                           rdflib.term.URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type'),  # noqa: E501
                                           rdflib.term.URIRef('http://purl.allotrope.org/ontologies/datacube-hdf-map#FunctionScaleMapping'))))  # noqa: E501
        # Step 3: Check the Index Function (contains type of function and parameters)
        # Step 4: Extract parameters
        if(mapping_function):
            param1 = float(list(g.objects(
                list(g.objects(mapping,
                               rdflib.term.URIRef('http://purl.allotrope.org/ontologies/datacube-hdf-map#indexFunction')))[0],  # noqa: E501
                rdflib.term.URIRef('http://purl.allotrope.org/ontologies/datacube-hdf-map#parameter1')))[0])  # noqa: E501
            param2 = float(list(g.objects(
                list(g.objects(mapping,
                               rdflib.term.URIRef('http://purl.allotrope.org/ontologies/datacube-hdf-map#indexFunction')))[0],  # noqa: E501
                rdflib.term.URIRef('http://purl.allotrope.org/ontologies/datacube-hdf-map#parameter2')))[0])  # noqa: E501
            return param1, param2
    return param1, param2


def read_adf_description(path, wl_len):
    """
    Queries the adf data description layer to extract the wavlength vector.
    For this query, the h5ld package is required which can be installed by
    editable pip install from https://github.com/laura-dirocco/h5ld. In case there
    are problems with installation, the user can give start and stop values
    as set on the DAD manually.
    """
    try:
        wl_slope, wl_start = get_function_paramenters(path)
        wl_stop = wl_start + wl_len * wl_slope
    except AttributeError:
        print("If the h5ld package cannot be installed on your machine, "
              "you have to give the wavelength values manually for adf data.")
    if wl_start == 0 and wl_slope == 0:
        wl_start = 190
        wl_stop = wl_start + wl_len
    wavelength = list(np.linspace(wl_start, wl_stop, wl_len))
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


def read_adf(path, wl_high_pass=None, wl_low_pass=None):
    """
    Reads adf files as exported by the Agilent ADF Adapter.
    """
    absorbance, time = read_adf_datacube(path)
    wavelength = read_adf_description(path, absorbance.shape[0])

    df = pd.DataFrame(np.swapaxes(absorbance, 0, 1), columns=wavelength)
    df.insert(0, "time", time)
    df = preprocess_df(df)

    df = pd.melt(df, id_vars='time', value_vars=df.columns[1:],
                 var_name='wavelength', value_name='absorbance')
    df['wavelength'] = df['wavelength'].astype(float)
    df = apply_filter(df, wl_high_pass, wl_low_pass)
    data, time, wavelength = df_to_array(df)
    return data, time, wavelength

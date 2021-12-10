#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 14:10:08 2021

@author: haascp
"""
import os

from mocca.dad_data.apis.chemstation_api import read_csv_agilent, tidy_df_agilent

def test_read_csv_agilent():
    file_path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(file_path,
                             "example_data/chemstation_api/test_data_chemstation")
    df = read_csv_agilent(data_path)
    assert len(df) == 6603
    assert len(df.columns) == 347

def test_tidy_df_agilent():
    file_path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(file_path,
                             "example_data/chemstation_api/test_data_chemstation")
    df = read_csv_agilent(data_path)
    df_tidy = tidy_df_agilent(df)
    assert len(df_tidy.columns) == 3

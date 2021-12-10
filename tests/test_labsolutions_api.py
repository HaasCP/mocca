#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 14:34:27 2021

@author: haascp
"""
import os

from mocca.dad_data.apis.labsolutions_api import read_txt_shimadzu

def test_read_txt_shimadzu():
    file_path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(file_path,
                             "example_data/labsolutions_api/labsolutions_test_data.txt")
    df = read_txt_shimadzu(data_path)
    assert len(df.columns) == 3
    
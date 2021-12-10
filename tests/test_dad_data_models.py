#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 13:59:54 2021

@author: haascp
"""
import os
import numpy as np

from mocca.dad_data.models import DadData, GradientData, CompoundData
from mocca.dad_data.utils import get_compound_names, get_compound_concentration

def test_dad_data_chemstation_api():
    file_path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(file_path,
                             "example_data/chemstation_api/test_data_chemstation")
    dad_data = DadData("chemstation", data_path)
    assert dad_data.path == data_path
    assert dad_data.hplc_system_tag == "chemstation"
    assert len(dad_data.time) == 6603
    assert len(dad_data.wavelength) == 342  # four less than raw data due to filter
    assert dad_data.data.shape == (342, 6603)
    assert dad_data.detector_limit == 2000

def test_dad_data_labsolutions_api():
    file_path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(file_path,
                             "example_data/labsolutions_api/labsolutions_test_data.txt")
    dad_data = DadData("labsolutions", data_path)
    assert dad_data.path == data_path
    assert dad_data.hplc_system_tag == "labsolutions"
    assert len(dad_data.time) == 3001
    assert len(dad_data.wavelength) == 160  # four less than raw data due to filter
    assert dad_data.data.shape == (160, 3001)
    assert dad_data.detector_limit == 2000

def test_gradient_data():
    file_path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(file_path,
                             "example_data/labsolutions_test_gradient.txt")
    dad_data = DadData("labsolutions", data_path)
    gradient_data = GradientData("labsolutions", data_path)
    assert gradient_data.path == data_path
    assert gradient_data.hplc_system_tag == "labsolutions"
    assert len(gradient_data.time) == 3001
    assert len(gradient_data.wavelength) == 160  # four less than raw data due to filter
    assert gradient_data.detector_limit == 2000
    assert gradient_data.data.shape == dad_data.data.shape
    assert not (gradient_data.data==dad_data.data).all()

def test_compound_data():
    file_path = os.path.dirname(os.path.abspath(__file__))
    gradient_path = os.path.join(file_path,
                                 "example_data/labsolutions_test_gradient.txt")
    gradient_data = GradientData("labsolutions", gradient_path)
    data_path = os.path.join(file_path,
                             "example_data/labsolutions_api/labsolutions_test_data.txt")
    test_input = {'a': 0, 'b': 0.1, 'c': 0, 'd': 1}
    compound_data = CompoundData("labsolutions", data_path, test_input, gradient_data)
    assert compound_data.path == data_path
    assert compound_data.hplc_system_tag == "labsolutions"
    assert len(compound_data.time) == 3001
    assert len(compound_data.wavelength) == 160  # four less than raw data due to filter
    assert compound_data.detector_limit == 2000
    assert np.amin(compound_data.data) > -20
    assert get_compound_names(compound_data) == ['b', 'd']
    assert get_compound_concentration(compound_data, 'a') == 0
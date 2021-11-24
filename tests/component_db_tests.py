#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 10:32:31 2021

@author: haascp
"""

import pytest


# COMPONENT DATABASE CLASS TESTS


def test_component_database_1():
    # try adding peaks to component database
    peak_1 = Peak(left=130, right=170, maximum=150, dataset=test_data_1)
    peak_2 = Peak(left=280, right=330, maximum=300, dataset=test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')
    component_database.add_peak(peak_2, 'Component 2')
    assert 'Component 1' in component_database
    assert 'Component 2' in component_database
    assert 'Component 3' not in component_database
    assert 'Component' not in component_database
    assert component_database['Component 1'].left == 130
    assert component_database['Component 1'].right == 170
    assert component_database['Component 1'].maximum == 150
    assert np.corrcoef(component_database['Component 1'].spectra,
                       test_spectra[1])[1, 0] > 0.999  # check that spectra are same
    components = []
    for component in component_database:
        components.append(component)
    assert len(components) == 2


def test_component_database_2():
    # component database add peak with same name
    with pytest.raises(Exception):
        peak_1 = Peak(left=130, right=150, maximum=100, dataset=test_data_1)
        peak_2 = Peak(left=280, right=330, maximum=300, dataset=test_data_1)
        component_database = ComponentDatabase()
        component_database.add_peak(peak_1, 'Component 1')
        component_database.add_peak(peak_2, 'Component 1')


def test_component_database_3():
    # component database get unadded item
    with pytest.raises(Exception):
        peak_1 = Peak(left=130, right=150, maximum=100, dataset=test_data_1)
        peak_2 = Peak(left=280, right=330, maximum=300, dataset=test_data_1)
        component_database = ComponentDatabase()
        component_database.add_peak(peak_1, 'Component 1')
        component_database.add_peak(peak_2, 'Component 2')
        component_database['Component 3']
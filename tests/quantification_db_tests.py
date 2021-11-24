#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 10:33:19 2021

@author: haascp
"""

import pytest

def test_peak_quantification_1():
    peak_1 = Peak(left=90, right=110, maximum=100, dataset=test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')

    qd = QuantificationDatabase()

    # error because self.integral and self.compound_id is None
    with pytest.raises(Exception):
        peak_1.quantify_peak(qd)

    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(1, 0.9))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(2, 2.1))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(3, 3))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(4, 4))

    peak_1.integrate_peak()
    print(peak_1.integral)
    peak_1.check_peak(detector_limit=10)
    peak_1.set_compound_id(component_database)
    peak_1.quantify_peak(qd)
    # quant database formula should be roughly y = x
    assert abs(peak_1.concentration - peak_1.integral) < 0.1


def test_peak_quantification_2():
    peak_1 = Peak(left=90, right=110, maximum=100, dataset=test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')

    qd = QuantificationDatabase()
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(1, 2))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(2, 5))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(3, 8))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(4, 10))

    peak_1.integrate_peak()
    peak_1.check_peak(detector_limit=10)
    peak_1.set_compound_id(component_database)
    peak_1.quantify_peak(qd)
    # quant database formula should be roughly y = 2.5x
    assert abs(peak_1.concentration - 2.5 * peak_1.integral) < 0.1


def test_peak_quantification_3():
    peak_1 = Peak(left=90, right=110, maximum=100, dataset=test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')

    qd = QuantificationDatabase()
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(1, 5))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(2, 5))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(3, 5))
    qd.add_compound_concentration(compound_id=('Component 1'), conc_info=(4, 5))

    peak_1.integrate_peak()
    peak_1.check_peak(detector_limit=10)
    peak_1.set_compound_id(component_database)
    peak_1.quantify_peak(qd)
    # make sure regression line passes through 0 still
    assert abs(peak_1.concentration - 5) > 0.1


def test_peak_quantification_4():
    peak_1 = Peak(left=90, right=110, maximum=100, dataset=test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')

    qd = QuantificationDatabase()
    qd.add_compound_concentration(compound_id=('Nonexistent 1'), conc_info=(1, 5))

    peak_1.integrate_peak()
    peak_1.check_peak(detector_limit=10)
    peak_1.set_compound_id(component_database)
    with pytest.raises(Exception):  # test nonexistent compound error
        peak_1.quantify_peak(qd)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 10:14:38 2021

@author: haascp
"""

# tests for peak_db.py
from mocca.peak_db import PeakDatabase
from mocca.peak import Peak
from chromatogram_gen import generate_test_chromatograms, plot_test_data

import logging
import pytest
# TODO: make dataset actually point to an object of DADData class,
# rather than Chromatogram class

# GET TEST DATA
test_data = generate_test_chromatograms()

# TEST DATA VISUALIZATION
"""
for data in test_data:
    plot_test_data(data)
"""

# CREATE LOGGER
LOGGER = logging.getLogger(__name__)

# ACTUAL TESTS
def test_init():
    test_db = PeakDatabase()
    assert isinstance(test_db.peaks, list)

def test_insert_1():  # insert new peak
    test_db = PeakDatabase()
    peak_1 = Peak(left=100, right=200, maximum=100, dataset=test_data[0])
    peak_2 = Peak(left=250, right=350, maximum=300, dataset=test_data[0])
    test_db.insert_peak(peak_1)
    test_db.insert_peak(peak_2)
    assert peak_1 in test_db.peaks
    assert peak_2 in test_db.peaks

def test_insert_2(caplog):
    test_db = PeakDatabase()
    peak_1 = Peak(left=100, right=200, maximum=100, dataset=test_data[0])
    peak_2 = Peak(left=100, right=200, maximum=100, dataset=test_data[0], saturation=True)
    with caplog.at_level(logging.WARNING):
        test_db.insert_peak(peak_1)
        test_db.insert_peak(peak_2)
    assert len(test_db.peaks)==1
    assert test_db.peaks[0].saturation==True
    assert "Warning: Peak" in caplog.text

def test_insert_3(caplog):
    test_db = PeakDatabase()
    peak_1 = Peak(left=100, right=200, maximum=100, dataset=test_data[0])
    peak_2 = {"left": 100, "right": 200, "maximum": 100, "dataset": test_data[0]}
    with caplog.at_level(logging.WARNING):
        test_db.insert_peak(peak_1)
        test_db.insert_peak(peak_2)
    assert len(test_db.peaks)==1
    assert "Warning: Data" in caplog.text

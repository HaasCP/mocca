#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 10:14:38 2021

@author: haascp
"""

# tests for peak_db.py
from mocca.databases.peak_db import PeakDatabase
from mocca.peak.models import ProcessedPeak, PickedPeak
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

def test_iter():
    test_db = PeakDatabase()
    peak_1 = ProcessedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    peak_2 = ProcessedPeak(left=250, right=350, maximum=300, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    test_db.peaks = [peak_1, peak_2]
    for i, peak in enumerate(test_db):
        assert peak == test_db.peaks[i]

def test_contains():
    test_db = PeakDatabase()
    peak_1 = ProcessedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    peak_2 = ProcessedPeak(left=250, right=350, maximum=300, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    test_db.peaks = [peak_1]
    assert peak_1 in test_db
    assert not peak_2 in test_db

def test_insert_1():  # insert new peak
    test_db = PeakDatabase()
    peak_1 = ProcessedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    peak_2 = ProcessedPeak(left=250, right=350, maximum=300, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    test_db.insert_peak(peak_1)
    test_db.insert_peak(peak_2)
    assert peak_1 in test_db.peaks
    assert peak_2 in test_db.peaks

def test_insert_2(caplog):
    test_db = PeakDatabase()
    peak_1 = ProcessedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    peak_2 = ProcessedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id="a", integral=12,
                            concentration=12.3)
    with caplog.at_level(logging.WARNING):
        test_db.insert_peak(peak_1)
        test_db.insert_peak(peak_2)
    assert len(test_db.peaks)==1
    assert test_db.peaks[0].compound_id=="a"
    assert "Warning: Peak" in caplog.text

def test_insert_3(caplog):
    test_db = PeakDatabase()
    peak_1 = ProcessedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    peak_2 = PickedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1)
    test_db.insert_peak(peak_1)
    with pytest.raises(Exception):
        test_db.insert_peak(peak_2)
    assert len(test_db.peaks)==1

#        logging.warning("{}".format(test_db.peaks))
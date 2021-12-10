#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 10:49:37 2021

@author: haascp
"""

from mocca.peak.models import ProcessedPeak
from mocca.components.databases import QualiComponentDatabase
from mocca.peak.database import PeakDatabase
from mocca.components.utils import get_filtered_peaks_by_compound, average_spectra_over_peaks

from chromatogram_gen import generate_test_chromatograms, plot_test_data

import logging
import pytest

from mocca.components.models import QualiComponent
from mocca.components.funcs import create_quali_component


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

def create_test_peak_db():
    test_db = PeakDatabase()
    db_peaks = []
    for i in range(11):
        peak = ProcessedPeak(left=100+i, right=200+i, maximum=100+i, dataset=test_data[0], idx=1,
                            saturation=False, pure=True, compound_id=str(int(i/3)), integral=12,
                            concentration=12.3)
        db_peaks.append(peak)
    test_db.peaks = db_peaks
    return test_db

def test_create_quali_component_1():
    test_db = create_test_peak_db()
    peak_dict = get_filtered_peaks_by_compound(test_db, filter_function=None)
    quali_components = []
    for key, value in peak_dict.items():
        comp = create_quali_component(value)
        quali_components.append(comp)
    for comp in quali_components:
        assert isinstance(comp, QualiComponent)
    assert quali_components[0].left == 101
    assert quali_components[1].spectrum == average_spectra_over_peaks(test_db.peaks[3:6])

def test_create_quali_component_2():
    test_db = create_test_peak_db()
    with pytest.raises(Exception):
         create_quali_component(test_db.peaks)

def create_test_components():
    test_db = create_test_peak_db()
    peak_dict = get_filtered_peaks_by_compound(test_db, filter_function=None)
    quali_components = []
    for key, value in peak_dict.items():
        comp = create_quali_component(value)
        quali_components.append(comp)
    return quali_components

def test_iter():
    quali_db = QualiComponentDatabase()
    quali_db.items = [0, 1, 2, 3, 4]
    for i, item in enumerate(quali_db):
        assert item == i

def test_contains():
    quali_db = QualiComponentDatabase()
    quali_db.items = create_test_components()
    assert "2" in quali_db

def test_getitem():
    quali_db = QualiComponentDatabase()
    quali_db.items = create_test_components()
    for comp in quali_db.items:
        if comp.compound_id == "2":
            test = comp
    assert quali_db["2"] == test

def test_insert_item_1():
    quali_db = QualiComponentDatabase()
    comps = create_test_components()
    quali_db.items = comps[:2]
    quali_db.insert_item(comps[2], unique=True)
    assert len(quali_db.items) == 3
    with pytest.raises(Exception):
        assert quali_db.insert_item(comps[2], unique=True)
    quali_db.insert_item(comps[2], unique=False)
    assert len(quali_db.items) == 4

def test_delete_all_items():
    quali_db = QualiComponentDatabase()
    quali_db.items = create_test_components()
    quali_db.delete_all_items()
    assert len(quali_db.items) == 0
    assert isinstance(quali_db.items, list)

def test_update():
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db)
    assert len(quali_db.items) == 4
    compound_ids = ["0", "1", "2", "3"]
    for compound_id in compound_ids:
        assert compound_id in quali_db

def filter_func(peaks):
    return [peak for peak in peaks if int(peak.compound_id) > 1]

def test_update_3():
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=filter_func)
    assert len(quali_db.items) == 2
    compound_ids = ["0", "1"]
    for compound_id in compound_ids:
        assert not compound_id in quali_db
        compound_ids = ["2", "3"]
        for compound_id in compound_ids:
            assert compound_id in quali_db

def test_insert_by_compound_id_1():
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=filter_func)
    quali_db.insert_by_compound_id(peak_db, "1", peak_filter_function=filter_func)
    assert len(quali_db.items) == 2
    compound_ids = ["0", "1"]
    for compound_id in compound_ids:
        assert not compound_id in quali_db
        compound_ids = ["2", "3"]
        for compound_id in compound_ids:
            assert compound_id in quali_db
    quali_db.insert_by_compound_id(peak_db, "1", peak_filter_function=None)
    assert len(quali_db.items) == 3
    compound_ids = ["0"]
    for compound_id in compound_ids:
        assert not compound_id in quali_db
        compound_ids = ["1", "2", "3"]
        for compound_id in compound_ids:
            assert compound_id in quali_db

def test_insert_by_compound_id_2():
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=filter_func)
    quali_db.insert_by_compound_id(peak_db, "1", peak_filter_function=None)
    quali_db.insert_by_compound_id(peak_db, "1", peak_filter_function=None)
    assert len(quali_db.items) == 3
    compound_ids = ["0"]
    for compound_id in compound_ids:
        assert not compound_id in quali_db
        compound_ids = ["1", "2", "3"]
        for compound_id in compound_ids:
            assert compound_id in quali_db
# TODO: Write tests after we can assign compound_id to peaks
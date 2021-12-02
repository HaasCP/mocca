#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 10:49:37 2021

@author: haascp
"""

# tests for peak_db.py
from mocca.component import Component, ComponentDatabase
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
    test_db = ComponentDatabase()
    assert isinstance(test_db.components, list)
    
# TODO: Write tests after we can assign compound_id to peaks
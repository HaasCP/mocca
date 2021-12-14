#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 17:57:52 2021

@author: haascp
"""

"""
from mocca.peak.preprocessor import preprocess_peak

def test_preprocess_peak():
    peak = PickedPeak(left=99, right=101, maximum=100, dataset=test_data[0], idx=1)
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    
    preprocessed_peak = preprocess_peak(peak, quali_db, absorbance_threshold=0.1, detector_limit=10,
    spectrum_correl_thresh= 0.99999, relative_distance_thresh=0.1)
    
    assert preprocessed_peak.left < 99 and preprocessed_peak.right > 101
    assert not preprocessed_peak.saturation
    assert preprocessed_peak.pure
    assert preprocessed_peak.integral > 0.9
    assert len(preprocessed_peak.matches) == 2
    assert preprocessed_peak.matches[0]['compound_id'] == "0" and preprocessed_peak.matches[1]['compound_id'] == "1"
"""
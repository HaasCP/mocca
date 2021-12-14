#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:15:37 2021

@author: haascp
"""
import logging

from chromatogram_gen import generate_test_chromatograms, plot_test_data
from mocca.components.databases import QualiComponentDatabase
from mocca.peak.database import PeakDatabase

# MODIFY AS NEEDED FOR TESTING
show_peak_purity_analytics = False

# GET TEST DATA
test_data = generate_test_chromatograms()

# TEST DATA VISUALIZATION
"""
for data in test_data:
    plot_test_data(data)
"""

# ACTUAL TESTS
from mocca.peak.models import (PickedPeak, CheckedPeak, IntegratedPeak, 
                               CorrectedPeak, PreprocessedPeak, ProcessedPeak)


# TEST PEAK PROCESS FUNCTIONS
from mocca.peak.expand import expand_peak

def test_peak_expand_1():
    # expand big peak
    peak = PickedPeak(left=848, right=852, maximum=850, dataset=test_data[0], idx=1)
    expanded_peak = expand_peak(peak, absorbance_threshold=0.1)
    left, right = (expanded_peak.left, expanded_peak.right)
    assert left < 835
    assert right > 865

def test_peak_expand_2():
    # expand small peak
    peak = PickedPeak(left=298, right=302, maximum=300, dataset=test_data[0], idx=1)
    expanded_peak = expand_peak(peak, absorbance_threshold=0.1)
    left, right = (expanded_peak.left, expanded_peak.right)
    assert left < 295
    assert right > 305

def test_peak_expand_3():
    # does not expand correct peak
    peak = PickedPeak(left=285, right=315, maximum=300, dataset=test_data[0], idx=1)
    expanded_peak = expand_peak(peak, absorbance_threshold=0.1)
    left, right = (expanded_peak.left, expanded_peak.right)
    assert abs(left - 285) <= 1  # +/- 1 error for tolerance due to noise
    assert abs(right - 315) <= 1

def test_peak_expand_4():
    # does not expand correct peak, no noise dataset
    peak = PickedPeak(left=285, right=315, maximum=300, dataset=test_data[6], idx=1)
    expanded_peak = expand_peak(peak, absorbance_threshold=0.1)
    left, right = (expanded_peak.left, expanded_peak.right)
    assert left == 285
    assert right == 315


from mocca.peak.check import check_peak_saturation, check_peak_purity, check_peak

def test_peak_saturation_1():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    assert not check_peak_saturation(peak, detector_limit=10)


def test_peak_saturation_2():
    peak = PickedPeak(left=150, right=250, maximum=200, dataset=test_data[2], idx=1)
    assert check_peak_saturation(peak, detector_limit=10)

def test_peak_purity_1():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1) 
    assert check_peak_purity(peak, show_analytics=show_peak_purity_analytics)

def test_peak_purity_2():
    # check for impure peak
    peak = PickedPeak(left=40, right=65, maximum=53, dataset=test_data[3], idx=1)
    assert not check_peak_purity(peak, show_analytics=show_peak_purity_analytics)

def test_peak_purity_3():
    # check for impure peak
    peak = PickedPeak(left=40, right=70, maximum=55, dataset=test_data[4], idx=1)
    assert not check_peak_purity(peak, show_analytics=show_peak_purity_analytics)

def test_peak_purity_4():
    # check for impure and saturated peak
    peak = PickedPeak(left=40, right=70, maximum=55, dataset=test_data[5], idx=1)
    assert not check_peak_purity(peak, show_analytics=show_peak_purity_analytics)

def test_check_peak():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    checked_peak = check_peak(peak, detector_limit=10, show_analytics=show_peak_purity_analytics)
    assert not checked_peak.saturation
    assert checked_peak.pure

from mocca.peak.integrate import integrate_peak

def test_integrate_peak_1():
    # check that peak integral is approximately 0 over empty peak area
    peak = CheckedPeak(left=550, right=650, maximum=600, dataset=test_data[0], idx=1,
                       pure=True, saturation=False)
    integrated_peak = integrate_peak(peak)
    assert integrated_peak.integral is not None
    assert integrated_peak.integral < 0.005  # approximate noise level

def test_integrate_peak_2():
    # check that peak integral is large over actual peak
    peak = CheckedPeak(left=90, right=110, maximum=100, dataset=test_data[0], idx=1,
                       pure=False, saturation=False)
    integrated_peak = integrate_peak(peak)
    assert integrated_peak.integral is not None
    assert integrated_peak.integral > 0.9
    # each peak should integrate to roughly its concentration in the chromatogram_gen
    # as spectra are normalized to area 1


from mocca.peak.correct import correct_offset

def test_correct_offset():
    peak = IntegratedPeak(left=90, right=110, maximum=100, dataset=test_data[0], idx=1,
                       pure=False, saturation=False, integral=0.99)
    corrected_peak_1 = correct_offset(peak, 5)
    corrected_peak_2 = correct_offset(peak, 0)
    corrected_peak_3 = correct_offset(peak, -10)
    assert corrected_peak_1.left == 85
    assert corrected_peak_2.right == 110
    assert corrected_peak_3.maximum == 110
    assert hasattr(corrected_peak_1, 'offset')

from mocca.peak.match import (get_spectrum_correl_coef, get_relative_distance,
                               get_similarity_dicts, get_filtered_similarity_dicts,
                               match_peak)

def create_test_peak_db():
    test_db = PeakDatabase()
    db_peaks = []
    for i in range(11):
        peak = ProcessedPeak(left=79+i, right=119+i, maximum=99+i, dataset=test_data[0], idx=1,
                            saturation=False, pure=True, compound_id=str(int(i/3)), integral=12, offset=0,
                            concentration=12.3)
        db_peaks.append(peak)
    test_db.peaks = db_peaks
    return test_db

def test_get_spectrum_correl_coef():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    correls = []
    for comp in quali_db:
        correls.append(get_spectrum_correl_coef(peak, comp))
    for i in range(len(correls)-1):
        assert correls[i] > correls[i+1]

def test_get_relative_distance():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    correls = []
    for comp in quali_db:
        correls.append(get_relative_distance(peak, comp))
    assert correls[0] == 0.0
    for i in range(len(correls)-1):
        assert correls[i] < correls[i+1]

def test_get_similarity_dict():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    dicts = get_similarity_dicts(peak, quali_db)
    for i in range(len(dicts)-1):
        assert dicts[i]['spectrum_correl_coef'] > dicts[i+1]['spectrum_correl_coef']
        assert dicts[i]['relative_distance'] < dicts[i+1]['relative_distance']
    assert dicts[0]['compound_id'] == "0"

def test_get_filtered_similarity_dicts_1():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    matches = get_filtered_similarity_dicts(peak, quali_db, spectrum_correl_coef_thresh=0.99999,
                                      relative_distance_thresh=0.0, print_out=True)
    assert len(matches) == 1
    assert matches[0]['compound_id'] == "0"

def test_get_filtered_similarity_dicts_2():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    matches = get_filtered_similarity_dicts(peak, quali_db, spectrum_correl_coef_thresh=0.99999,
                                      relative_distance_thresh=0.1, print_out=True)
    assert len(matches) == 2
    assert matches[0]['compound_id'] == "0" and matches[1]['compound_id'] == "1"

def test_match_peak_1():
    peak = CorrectedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1,
                          pure=True, saturation=True, integral=123, offset=0)
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    matched_peak = match_peak(peak, quali_db, spectrum_correl_coef_thresh=0.99999,
                              relative_distance_thresh=0.1, print_similarity_dicts=True)
    assert len(matched_peak.matches) == 2
    assert matched_peak.matches[0]['compound_id'] == "0" and matched_peak.matches[1]['compound_id'] == "1"

def test_match_peak_2():
    peak = CorrectedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1,
                          pure=False, saturation=True, integral=123, offset=0)
    peak_db = create_test_peak_db()
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    matched_peak = match_peak(peak, quali_db, spectrum_correl_coef_thresh=0.99999,
                              relative_distance_thresh=0.1, print_similarity_dicts=True)
    assert matched_peak.matches is None

    #logging.warning("{}".format(dicts))


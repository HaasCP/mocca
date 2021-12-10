#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 09:21:06 2021

@author: haascp
"""
from mocca.peak.models import PickedPeak, PreprocessedPeak
from mocca.components.databases import QualiComponentDatabase

from mocca.peak.expand import expand_peak
from mocca.peak.check import check_peak
from mocca.peak.integrate import integrate_peak
from mocca.peak.match import match_peak


def preprocess_peak(picked_peak: PickedPeak, component_db: QualiComponentDatabase,
                    absorbance_threshold: int, detector_limit: int,
                    spectrum_correl_thresh: float = 0.9,  # spectrum_correl_coef_thresh
                    relative_distance_thresh: float = 0.01,  # relative_distance_thresh
                    print_purity_check: bool = False,
                    print_compound_prediction: bool = False) -> PreprocessedPeak:
    """
    Peak preprocessing routine. Takes a picked peak and returns a preprocessed
    peak which was expanded, checked for saturation and purity, integrated
    and assigend with possible component matches.
    """
    expanded_peak = expand_peak(picked_peak, absorbance_threshold)
    checked_peak = check_peak(expanded_peak, detector_limit,
                              show_analytics=print_purity_check)
    integrated_peak = integrate_peak(checked_peak)
    preprocessed_peak = match_peak(integrated_peak, component_db,
                                   spectrum_correl_thresh,
                                   relative_distance_thresh,
                                   print_compound_prediction)
    return preprocessed_peak

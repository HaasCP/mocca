#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 16:08:40 2021

@author: haascp
"""

from mocca.peak.expand import expand_peak
from mocca.peak.check import check_peak
from mocca.peak.integrate import integrate_peak
from mocca.chromatogram.correct import correct_istd_offset
from mocca.decomposition.utils import check_any_compound_overlap
from mocca.decomposition.utils import check_absorbance_thresh
from mocca.decomposition.iterative_parafac import iterative_parafac
from mocca.peak.match import match_peak


def preprocess_chromatogram(chromatogram, quali_comp_db,
                            absorbance_threshold, detector_limit,
                            spectrum_correl_thresh, relative_distance_thresh,
                            print_purity_check=True,
                            print_compound_prediction=False,
                            print_parafac_analytics=False):
    """
    Preprocesses the chromatogram of picked peaks. It includes expanding,
    checking, integrating, correcting, resolving impures, and matching of the
    peaks in the chromatogram.
    """
    # 1. expand, 2. check, 3. integrate
    integrated_peaks = []
    for picked_peak in chromatogram.peaks:
        expanded_peak = expand_peak(picked_peak, absorbance_threshold)
        checked_peak = check_peak(expanded_peak, detector_limit,
                                  show_analytics=print_purity_check)
        integrated_peak = integrate_peak(checked_peak)
        integrated_peaks.append(integrated_peak)

    chromatogram.peaks = integrated_peaks

    # 4. correct
    chromatogram = correct_istd_offset(chromatogram, quali_comp_db,
                                       absorbance_threshold,
                                       spectrum_correl_thresh,
                                       relative_distance_thresh)
    # 5. resolve impure
    impure_peaks = [peak for peak in chromatogram if not peak.pure]
    relevant_impure_peaks = [peak for peak in impure_peaks if
                             check_any_compound_overlap(peak, quali_comp_db)]

    for impure_peak in relevant_impure_peaks:
        parafac_model = iterative_parafac(impure_peak, quali_comp_db,
                                          absorbance_threshold,
                                          relative_distance_thresh,
                                          spectrum_correl_thresh,
                                          print_parafac_analytics)

        chromatogram.parafac_models.append(parafac_model)

        for parafac_peak in parafac_model.peaks:
            if check_absorbance_thresh(parafac_peak, absorbance_threshold):
                chromatogram.peaks.append(parafac_peak)

    # 6. match
    matched_peaks = []
    for resolved_peak in chromatogram:
        new_peak = match_peak(resolved_peak, quali_comp_db,
                              spectrum_correl_thresh,
                              relative_distance_thresh,
                              print_compound_prediction)
        matched_peaks.append(new_peak)
    chromatogram.peaks = matched_peaks
    return chromatogram

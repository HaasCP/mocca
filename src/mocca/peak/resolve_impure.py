#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 08:59:36 2022

@author: haascp
"""
import numpy as np
import itertools

from mocca.peak.models import CorrectedPeak, IntegratedPeak
from mocca.dad_data.models import ParafacData
from mocca.dad_data.utils import sum_absorbance_by_time
from mocca.decomposition.iterative_parafac import iterative_parafac


def check_same_uvvis(parafac_factors, n_comps, spectrum_correl_coef_thresh):
    """
    Checks if any two parafac components share the same UV-Vis trace.
    """
    if all(np.corrcoef(parafac_factors[0][:, i], parafac_factors[0][:, j])[0, 1] > 
           spectrum_correl_coef_thresh for i, j in 
           itertools.combinations(list(range(n_comps)), 2)):
        return True
    else:
        return False


def create_pure_peak(impure_peak):
    """
    Takes an impure peak and returns its copy with the pure attribute True.
    """
    if type(impure_peak) == CorrectedPeak:
        parafac_peak = CorrectedPeak(left=impure_peak.left,
                                     right=impure_peak.right,
                                     maximum=impure_peak.maximum,
                                     offset=impure_peak.offset,
                                     dataset=impure_peak.dataset,
                                     idx=-impure_peak.idx,
                                     saturation=impure_peak.saturation,
                                     pure=True,
                                     integral=impure_peak.integral,
                                     istd=impure_peak.istd)
    elif type(impure_peak) == IntegratedPeak:
        parafac_peak = IntegratedPeak(left=impure_peak.left,
                                     right=impure_peak.right,
                                     maximum=impure_peak.maximum,
                                     offset=impure_peak.offset,
                                     dataset=impure_peak.dataset,
                                     idx=-impure_peak.idx,
                                     saturation=impure_peak.saturation,
                                     pure=True,
                                     integral=impure_peak.integral)
    else:
        raise TypeError(f"Given impure peak is of type {type(impure_peak)}. "
                        "Only mocca IntegratedPeak and CorrectedPeak types "
                        "are allowed.")
    return parafac_peak


def create_parafac_peak(impure_peak, parafac_comp_factors, boundaries, shift,
                        y_offset):
    """
    Return synthetic PARAFAC peaks created from the PARAFAC decomposition results.
    """
    # integral correction for y-shifted parafac data
    integral_correction = y_offset * len(impure_peak.dataset.wavelength) *\
        (boundaries[1] - boundaries[0])
    if type(impure_peak) == CorrectedPeak:
        parafac_peak = CorrectedPeak(left=boundaries[0],
                                     right=boundaries[1],
                                     maximum=(boundaries[0] +
                                              np.argmax(parafac_comp_factors[1])),
                                     offset=impure_peak.offset,
                                     dataset=ParafacData(impure_peak,
                                                         parafac_comp_factors,
                                                         boundaries,
                                                         shift,
                                                         y_offset),
                                     idx=-impure_peak.idx,
                                     saturation=impure_peak.saturation,
                                     pure=True,
                                     # reaction run in last slice of data tensor
                                     integral=parafac_comp_factors[2][-1] +\
                                         integral_correction,
                                     istd=impure_peak.istd)
    elif type(impure_peak) == IntegratedPeak:
        parafac_peak = IntegratedPeak(left=boundaries[0],
                                     right=boundaries[1],
                                     maximum=(boundaries[0] +
                                              np.argmax(parafac_comp_factors[1]) -
                                              shift),
                                     offset=impure_peak.offset,
                                     dataset=ParafacData(impure_peak,
                                                         parafac_comp_factors,
                                                         boundaries,
                                                         shift,
                                                         y_offset),
                                     idx=-impure_peak.idx,
                                     saturation=impure_peak.saturation,
                                     pure=True,
                                     # reaction run in last slice of data tensor
                                     integral=parafac_comp_factors[2][-1] +\
                                         integral_correction)
    else:
        raise TypeError(f"Given impure peak is of type {type(impure_peak)}. "
                        "Only mocca IntegratedPeak and CorrectedPeak types "
                        "are allowed.")
    return parafac_peak


def check_absorbance_thresh(parafac_peak, absorbance_threshold):
    """
    Checks if maximum absorbance in synthetically created PARAFAC peak dataset
    exceeds absorbance threshold.
    """
    max_absorbance = np.max(sum_absorbance_by_time(parafac_peak.dataset.data))
    return max_absorbance > absorbance_threshold


def create_parafac_peaks(impure_peak, parafac_factors, boundaries, iter_offset,
                         y_offset, absorbance_threshold,
                         spectrum_correl_coef_thresh=0.999):
    """
    Takes an impure peak and its PARAFAC decomposition results and creates
    new peaks in a readible form for the program. If two UV-Vis traces in the
    PARAFAC components are too similar, the peak is set to pure. Else, new
    PARAFAC peaks are created. Its dataset is created synthetically by using
    Spectra, Elution, and Integral generated from PARAFAC and filling the rest
    of the array up with zeros. PARAFAC peaks get an index of -impure_peak.idx
    """
    n_comps = parafac_factors[0].shape[1]

    if check_same_uvvis(parafac_factors, n_comps, spectrum_correl_coef_thresh):
        parafac_peak = create_pure_peak(impure_peak)
        return [parafac_peak], None

    else:
        if iter_offset > 0:
            shift = iter_offset
        else:
            shift = 0
        parafac_peaks = []
        chrom_peaks = []
        for i in range(n_comps):
            #  get factors for one parafac comonent
            parafac_comp_factors = (parafac_factors[0][:, i],
                                   parafac_factors[1][:, i],
                                   parafac_factors[2][:, i])
            parafac_peak = create_parafac_peak(impure_peak, parafac_comp_factors,
                                               boundaries, shift, y_offset)
            parafac_peaks.append(parafac_peak)
            
            if check_absorbance_thresh(parafac_peak, absorbance_threshold):
                chrom_peaks.append(parafac_peak)

        return chrom_peaks, (impure_peak, parafac_peaks)


def get_parafac_peaks(impure_peak, quali_comp_db, absorbance_threshold,
                      spectrum_correl_coef_thresh, relative_distance_thresh,
                      show_parafac_analytics):
    """
    Runs PARAFAC decomposition function and returns the calculated peaks as well
    as information required for the PARAFAC report.
    """
    parafac_factors, boundaries, iter_offset, y_offset =\
        iterative_parafac(impure_peak, quali_comp_db, relative_distance_thresh,
                          show_parafac_analytics)

    chrom_peaks, parafac_report_tuple =\
        create_parafac_peaks(impure_peak, parafac_factors, boundaries,
                             iter_offset, y_offset, absorbance_threshold,
                             spectrum_correl_coef_thresh)
    return chrom_peaks, parafac_report_tuple

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


def get_parafac_data_shift(iter_offset):
    if iter_offset > 0:
        shift = iter_offset
    else:
        shift = 0
    return shift


def create_parafac_peak(comp_i, parafac_model, iter_offset):
    """
    Return synthetic PARAFAC peaks created from the PARAFAC decomposition results.
    """
    shift = get_parafac_data_shift(iter_offset)
    parafac_comp_factors = (parafac_model.factors[0][:, comp_i],
                            parafac_model.factors[1][:, comp_i],
                            parafac_model.factors[2][:, comp_i])
    tensor = parafac_model.data_tensor
    impure_peak = parafac_model.impure_peak
    left_bound = tensor.boundaries[0]
    right_bound = tensor.boundaries[1]
    n_wls = len(parafac_model.impure_peak.dataset.wavelength)
    # integral correction for y-shifted parafac data
    integral_correction = tensor.y_offset * n_wls * (right_bound - left_bound)

    if type(impure_peak) == CorrectedPeak:
        parafac_peak = CorrectedPeak(left=left_bound,
                                     right=right_bound,
                                     maximum=(left_bound +
                                              np.argmax(parafac_comp_factors[1])),
                                     offset=impure_peak.offset,
                                     dataset=ParafacData(impure_peak,
                                                         parafac_comp_factors,
                                                         tensor.boundaries,
                                                         shift,
                                                         tensor.y_offset),
                                     idx=-impure_peak.idx,
                                     saturation=impure_peak.saturation,
                                     pure=True,
                                     # reaction run in last slice of data tensor
                                     integral=parafac_comp_factors[2][-1] +\
                                         integral_correction,
                                     istd=impure_peak.istd)
    elif type(impure_peak) == IntegratedPeak:
        parafac_peak = IntegratedPeak(left=left_bound,
                                     right=right_bound,
                                     maximum=(left_bound +
                                              np.argmax(parafac_comp_factors[1]) -
                                              shift),
                                     offset=impure_peak.offset,
                                     dataset=ParafacData(impure_peak,
                                                         parafac_comp_factors,
                                                         tensor.boundaries,
                                                         shift,
                                                         tensor.y_offset),
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


def check_same_uvvis(parafac_model, spectrum_correl_coef_thresh):
    """
    Checks if any two parafac components share the same UV-Vis trace.
    """
    spectra = parafac_model.factors[0]
    if all(np.corrcoef(spectra[:, i], spectra[:, j])[0, 1] > 
           spectrum_correl_coef_thresh for i, j in 
           itertools.combinations(list(range(parafac_model.n_comps)), 2)):
        return True
    else:
        return False


def create_parafac_peaks(parafac_model, iter_offset, absorbance_threshold,
                         spectrum_correl_coef_thresh):
    """
    Takes an impure peak and its PARAFAC decomposition results and creates
    new peaks in a readible form for the program. If two UV-Vis traces in the
    PARAFAC components are too similar, the peak is set to pure. Else, new
    PARAFAC peaks are created. Its dataset is created synthetically by using
    Spectra, Elution, and Integral generated from PARAFAC and filling the rest
    of the array up with zeros. PARAFAC peaks get an index of -impure_peak.idx
    """

    if check_same_uvvis(parafac_model, spectrum_correl_coef_thresh):
        parafac_peak = create_pure_peak(parafac_model.impure_peak)
        return [parafac_peak], None

    else:
        parafac_peaks = []
        chrom_peaks = []
        for i in range(parafac_model.n_comps):
            #  get factors for one parafac comonent
            parafac_peak = create_parafac_peak(i, parafac_model, iter_offset)
            parafac_peaks.append(parafac_peak)

            if check_absorbance_thresh(parafac_peak, absorbance_threshold):
                chrom_peaks.append(parafac_peak)

        return chrom_peaks, parafac_model


def get_parafac_peaks(impure_peak, quali_comp_db, absorbance_threshold,
                      spectrum_correl_coef_thresh, relative_distance_thresh,
                      show_parafac_analytics):
    """
    Runs PARAFAC decomposition function and returns the calculated peaks as well
    as information required for the PARAFAC report.
    """
    parafac_model, iter_offset, iter_objective_func =\
        iterative_parafac(impure_peak, quali_comp_db, relative_distance_thresh,
                          show_parafac_analytics)

    chrom_peaks, parafac_model =\
        create_parafac_peaks(parafac_model, iter_offset, absorbance_threshold,
                             spectrum_correl_coef_thresh)
    return chrom_peaks, parafac_model, iter_offset, iter_objective_func

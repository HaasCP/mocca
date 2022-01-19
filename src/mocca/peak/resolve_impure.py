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


def check_same_uvvis(parafac_factors, spectrum_correl_coef_thresh, n_comps):
    if any(np.corrcoef(parafac_factors[0][:, i], parafac_factors[0][:, j])[0, 1] > 
           spectrum_correl_coef_thresh for i, j in 
           itertools.combinations(list(range(n_comps)), 2)):
        return True
    else:
        return False

def create_parafac_peaks(impure_peak, parafac_factors, boundaries,
                         absorbance_threshold, spectrum_correl_coef_thresh):
    """
    Makes a new ProcessedPeak corresponding to PARAFAC-processed impure processed peak
    Spectra, Elution, and Integral generated from PARAFAC

    If the old peak had index i, then the new peak has index i_{component_num}
    """
    n_comps = parafac_factors[0].shape[1]

    if check_same_uvvis(parafac_factors, spectrum_correl_coef_thresh, n_comps):
        if type(impure_peak) == CorrectedPeak:
            parafac_peak = CorrectedPeak(left=impure_peak.left,
                                         right=impure_peak.right,
                                         maximum=impure_peak.maximum,
                                         dataset=impure_peak.dataset,
                                         idx=-impure_peak.idx,
                                         saturation=impure_peak.saturation,
                                         pure=True,
                                         integral=impure_peak.integral, # reaction run is last in run dimension
                                         offset=impure_peak.offset,
                                         istd=impure_peak.istd)
        elif type(impure_peak) == IntegratedPeak:
            parafac_peak = IntegratedPeak(left=impure_peak.left,
                                         right=impure_peak.right,
                                         maximum=impure_peak.maximum,
                                         dataset=impure_peak.dataset,
                                         idx=-impure_peak.idx,
                                         saturation=impure_peak.saturation,
                                         pure=True,
                                         integral=impure_peak.integral)
        else:
            raise TypeError(f"Given impure peak is of type {type(impure_peak)}. "
                            "Only mocca IntegratedPeak and CorrectedPeak types "
                            "are allowed.")
        return [parafac_peak], None

    else:
        parafac_peaks = []
        chrom_peaks = []
        integral_sum_comp_runs = []
        for i in range(n_comps):
            #  get factors for one parafac comonent
            parafac_comp_factors = (parafac_factors[0][:, i],
                                   parafac_factors[1][:, i],
                                   parafac_factors[2][:, i])
            if type(impure_peak) == CorrectedPeak:
                parafac_peak = CorrectedPeak(left=boundaries[0],
                                             right=boundaries[1],
                                             maximum=(boundaries[0] +
                                                      np.argmax(parafac_comp_factors[1])),
                                             dataset=ParafacData(impure_peak,
                                                                 parafac_comp_factors,
                                                                 boundaries),
                                             idx=-impure_peak.idx,
                                             saturation=impure_peak.saturation,
                                             pure=True,
                                             integral=parafac_comp_factors[2][-1], # reaction run is last in run dimension
                                             offset=0,
                                             istd=impure_peak.istd)
            elif type(impure_peak) == IntegratedPeak:
                parafac_peak = IntegratedPeak(left=boundaries[0],
                                             right=boundaries[1],
                                             maximum=(boundaries[0] +
                                                      np.argmax(parafac_comp_factors[1])),
                                             dataset=ParafacData(impure_peak,
                                                                 parafac_comp_factors,
                                                                 boundaries),
                                             idx=-impure_peak.idx,
                                             saturation=impure_peak.saturation,
                                             pure=True,
                                             integral=parafac_comp_factors[2][-1])
            else:
                raise TypeError(f"Given impure peak is of type {type(impure_peak)}. "
                                "Only mocca IntegratedPeak and CorrectedPeak types "
                                "are allowed.")
            parafac_peaks.append(parafac_peak)
            max_absorbance = np.max(sum_absorbance_by_time(parafac_peak.dataset.data))
            if max_absorbance > absorbance_threshold:
                chrom_peaks.append(parafac_peak)
            integral_sum_comp_runs.append(sum(parafac_comp_factors[2][:-1]))
    
        return chrom_peaks, (impure_peak, parafac_peaks)


def get_parafac_peaks(impure_peak, quali_comp_db, absorbance_threshold,
                      spectrum_correl_coef_thresh, show_parafac_analytics):
    """
    Runs PARAFAC decomposition function and returns the calculated peaks as well
    as information required for the PARAFAC report.
    """
    parafac_factors, boundaries = iterative_parafac(impure_peak, quali_comp_db)
    chrom_peaks, parafac_report_tuple = create_parafac_peaks(impure_peak,
                                                               parafac_factors,
                                                               boundaries,
                                                               absorbance_threshold,
                                                               spectrum_correl_coef_thresh)
    return chrom_peaks, parafac_report_tuple

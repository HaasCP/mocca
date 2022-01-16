#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 08:59:36 2022

@author: haascp
"""
import numpy as np

from mocca.peak.models import CorrectedPeak, IntegratedPeak
from mocca.dad_data.models import ParafacData
from mocca.dad_data.utils import sum_absorbance_by_time
from mocca.decomposition.iterative_parafac import iterative_parafac


def create_parafac_peaks(impure_peak, parafac_factors, boundaries,
                         absorbance_threshold):
    """
    Makes a new ProcessedPeak corresponding to PARAFAC-processed impure processed peak
    Spectra, Elution, and Integral generated from PARAFAC

    If the old peak had index i, then the new peak has index i_{component_num}
    """
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
        if np.max(sum_absorbance_by_time(parafac_peak.dataset.data)) > absorbance_threshold:
            chrom_peaks.append(parafac_peak)
        
        integral_sum_comp_runs.append(sum(parafac_comp_factors[2][:-1]))
        print(f"Summed integral for component {i}: {sum(parafac_comp_factors[2][:])}")

    #summed_data = get_peak_data(parafac_peaks[0])
    #if len(parafac_peaks) > 1:
    #    for parafac_peak in parafac_peaks[1:]:
    #        summed_data = np.add(summed_data, get_peak_data(parafac_peak))

    #print(np.argmax(np.sum(summed_data, axis=0)))
    #display(contour_map(summed_data, list(range(summed_data.shape[1])),
    #                    list(range(summed_data.shape[0]))))
    
    #impure_ze_data = get_impure_ze_peak(impure_peak, boundaries)
    
    #cov_matrix = np.cov(impure_ze_data, summed_data)
    #egnvalues, egnvectors = eigh(cov_matrix)
    #
    # Determine explained variance
    #
    #total_egnvalues = sum(egnvalues)
    #var_exp = [(i/total_egnvalues) for i in sorted(egnvalues, reverse=True)]
    #print(var_exp[0])
    
    return chrom_peaks, (impure_peak, parafac_peaks)


def get_parafac_peaks(impure_peak, quali_comp_db, absorbance_threshold,
                      show_parafac_analytics):
    parafac_factors, boundaries = iterative_parafac(impure_peak, quali_comp_db)
    chrom_peaks, parafac_report_tuple = create_parafac_peaks(impure_peak,
                                                               parafac_factors,
                                                               boundaries,
                                                               absorbance_threshold)
    return chrom_peaks, parafac_report_tuple
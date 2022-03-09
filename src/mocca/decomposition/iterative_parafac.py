#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 11:15:19 2022

@author: haascp
"""
from mocca.decomposition.parafac_funcs import parafac


def run_parafac_iter(offset, impure_peak, quali_comp_db, show_parafac_analytics):
    """
    One iteration of the PARAFAC routine plus returning the result of the
    objective function.
    """
    parafac_model = parafac(impure_peak, quali_comp_db, offset,
                            show_parafac_analytics)

    #impure_sum = get_impure_integral_sum(parafac_model.factors,
    #                                     show_parafac_analytics)

    #return impure_sum, parafac_model


def iterative_parafac(impure_peak, quali_comp_db, absorbance_threshold,
                      relative_distance_thresh, spectrum_correl_coef_thresh,
                      show_parafac_analytics):
    """
    The trilinearity-breaking mode retention time requires iterative PARAFAC
    algorithm.
    """
    # set peak offset as middle point of the iterator
    offset_seed = -impure_peak.offset

    len_iterator = int(relative_distance_thresh * len(impure_peak.dataset.time))
    offset_iterator = [i + offset_seed - len_iterator for i in
                       list(range(len_iterator * 2 + 1))]

    mse_opt = float("inf")

    iter_objective_func = []
    for offset in offset_iterator:
        parafac_model = parafac(impure_peak, quali_comp_db, offset,
                                show_parafac_analytics)
        iter_objective_func.append((offset, parafac_model.impure_mse))
        if parafac_model.impure_mse <= mse_opt:
            mse_opt = parafac_model.impure_mse
            parafac_model_opt = parafac_model

    parafac_model_opt.iter_objective_func = iter_objective_func
    parafac_model_opt.create_parafac_peaks(absorbance_threshold,
                                           spectrum_correl_coef_thresh)

    return parafac_model_opt

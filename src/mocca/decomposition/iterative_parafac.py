#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 11:15:19 2022

@author: haascp
"""
from mocca.decomposition.parafac_funcs import parafac


def get_impure_integral_sum(parafac_factors, show_parafac_analytics):
    """
    This objective function assumes that the PARAFAC model is best, when the
    summed integral of all components in the impure peak slice is maximized.
    """
    integrals = parafac_factors[2][-1, :]
    sum_i = sum(integrals)
    if show_parafac_analytics:
        print(f"impure_peak_sum = {sum_i}")
    return sum_i


def run_parafac_iter(offset, impure_peak, quali_comp_db, show_parafac_analytics):
    """
    One iteration of the PARAFAC routine plus returning the result of the
    objective function.
    """
    parafac_model = parafac(impure_peak, quali_comp_db, offset,
                            show_parafac_analytics)

    impure_sum = get_impure_integral_sum(parafac_model.factors,
                                         show_parafac_analytics)

    return impure_sum, parafac_model


def iterative_parafac(impure_peak, quali_comp_db, relative_distance_thresh,
                      spectrum_correl_coef_thresh, show_parafac_analytics):
    """
    The trilinearity-breaking mode retention time requires iterative PARAFAC
    algorithm.
    """
    # set peak offset as middle point of the iterator
    offset_seed = -impure_peak.offset

    len_iterator = int(relative_distance_thresh * len(impure_peak.dataset.time))
    offset_iterator = [i + offset_seed - len_iterator for i in
                       list(range(len_iterator * 2 + 1))]

    impure_sum_opt = 0

    iter_objective_func = []
    for offset in offset_iterator:
        impure_sum_new, parafac_model = run_parafac_iter(offset, impure_peak,
                                                         quali_comp_db,
                                                         show_parafac_analytics)
        iter_objective_func.append((offset, impure_sum_new))
        if impure_sum_new >= impure_sum_opt:
            impure_sum_opt = impure_sum_new
            parafac_model_opt = parafac_model

    parafac_model_opt.iter_objective_func = iter_objective_func
    parafac_model_opt.create_parafac_peaks(spectrum_correl_coef_thresh)

    return parafac_model_opt

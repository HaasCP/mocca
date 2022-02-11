#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 11:15:19 2022

@author: haascp
"""
import math
from mocca.decomposition.parafac_funcs import parafac


def get_non_comp_sum(parafac_factors, start_slice, end_slice):
    """
    Get summed integrals from all components except the highest inegral component.
    """
    n_comps = parafac_factors[0].shape[1]
    non_comp_sum_list = []
    for i in range(n_comps):
        sorted_integrals = sorted(parafac_factors[2][:, i][start_slice:end_slice])
        non_comp_sum_i = sum(sorted_integrals)
        non_comp_sum_list.append(non_comp_sum_i)
    sorted_sums_ex_max = sorted(non_comp_sum_list)[:-1]
    return sum(sorted_sums_ex_max)


def get_all_non_comp_sum(parafac_factors, comp_tensor_shape,
                         show_parafac_analytics):
    """
    Iterative PARAFAC approach is a minimization problem. Based on the component
    tensor data shape, the maximal sum of integrals is found for every component
    part in the data tensor. The integrals of all other components in this
    component part are summed up and should be ideally zero (only component
    spectrum contributes to component part of tensor). This final sum
    is the minimization objective for the iterative PARAFAC approach.
    """
    start_slice = 0
    end_slice = 0
    non_comp_sum_list = []
    for i in range(len(comp_tensor_shape)):
        end_slice += comp_tensor_shape[i]
        non_comp_sum = get_non_comp_sum(parafac_factors, start_slice, end_slice)
        non_comp_sum_list.append(non_comp_sum)
        start_slice += comp_tensor_shape[i]
    if show_parafac_analytics:
        print(f"objective_val_non_comp = {sum(non_comp_sum_list)}")
    return sum(non_comp_sum_list)


def get_comp_sum(parafac_factors, start_slice, end_slice):
    """
    Get summed integrals from the highest inegral component.
    """
    n_comps = parafac_factors[0].shape[1]
    non_comp_sum_list = []
    for i in range(n_comps):
        sorted_integrals = sorted(parafac_factors[2][:, i][start_slice:end_slice])
        non_comp_sum_i = sum(sorted_integrals)
        non_comp_sum_list.append(non_comp_sum_i)
    max_sum = sorted(non_comp_sum_list)[-1]
    return max_sum


def get_all_comp_sum(parafac_factors, comp_tensor_shape, show_parafac_analytics):
    start_slice = 0
    end_slice = 0
    comp_sum_list = []
    for i in range(len(comp_tensor_shape)):
        end_slice += comp_tensor_shape[i]
        comp_sum = get_comp_sum(parafac_factors, start_slice, end_slice)
        comp_sum_list.append(comp_sum)
        start_slice += comp_tensor_shape[i]
    if show_parafac_analytics:
        print(f"objective_val_comp = {sum(comp_sum_list)}")
    return sum(comp_sum_list)


def get_total_integral_sum(parafac_factors, show_parafac_analytics):
    n_comps = parafac_factors[0].shape[1]
    sum_list = []
    for i in range(n_comps):
        integrals = parafac_factors[2][:, i]
        sum_i = sum(integrals)
        sum_list.append(sum_i)
    total_sum = sum(sum_list)
    if show_parafac_analytics:
        print(f"objective_total_sum = {total_sum}")
    return total_sum


def get_impure_integral_sum(parafac_factors, show_parafac_analytics):
    integrals = parafac_factors[2][-1, :]
    sum_i = sum(integrals)
    if show_parafac_analytics:
        print(f"impure_peak_sum = {sum_i}")
    return sum_i


def offset_opt_func(iter_offset_new, impure_peak, quali_comp_db,
                    show_parafac_analytics):
    parafac_factors_new, boundaries_new, *_ = parafac(impure_peak,
                                                      quali_comp_db,
                                                      iter_offset_new,
                                                      show_parafac_analytics)

    impure_sum_new = get_impure_integral_sum(parafac_factors_new,
                                         show_parafac_analytics)
    return impure_sum_new, parafac_factors_new, boundaries_new


def iterative_parafac(impure_peak, quali_comp_db, relative_distance_thresh,
                      show_parafac_analytics):
    """
    The trilinearity-breaking mode retention time requires iterative PARAFAC
    algorithm. The sum of the non-component integrals in the component part
    of the data tensor is the minimization objective upon changing the offset
    in the retention time mode. Empirically it was found that this value has
    two maxima and one local minimum around the true offset with the minimum
    exactly at the correct offset. If the offset gets very large, the value
    gets minimized again. The algorithm assumes the retention time dimension
    to be reproducible enough that the start point with offset == 0 lies somewhere
    between the two maxima. It moves first in one direction. If the slope is
    positive it reverses direction. At the point where slope is positive in
    both directions, the best offset is reached.
    """
    # start point is offset to have good initial guess
    offset_init = -impure_peak.offset
    # comp_tensor_shape should not change over iterative parafac runs
    # initialize optimization problem
    parafac_factors_cur, boundaries_cur, comp_tensor_shape, y_offset =\
        parafac(impure_peak, quali_comp_db, offset_init, show_parafac_analytics)

    offset_opt = offset_init
    impure_sum_opt = get_impure_integral_sum(parafac_factors_cur,
                                             show_parafac_analytics)
    
    neg_term_count = 0  # counts not improving iterations in negative dir
    pos_term_count = 0  # counts not improving iterations in positive dir
    offset = 0  # determines step size away from offset_init in both dirs
    termination_criterion = False  # optimized when True
    
    while not termination_criterion:        
        offset += 1
        offset_neg = offset_init - offset
        offset_pos = offset_init + offset
        if neg_term_count < 3:
            impure_sum_new, parafac_factors_new, boundaries_new =\
                offset_opt_func(offset_neg, impure_peak, quali_comp_db,
                                show_parafac_analytics)
            if impure_sum_new > impure_sum_opt:
                impure_sum_opt = impure_sum_new
                offset_opt = offset_neg
                parafac_factors_cur = parafac_factors_new
                boundaries_cur = boundaries_new
                neg_term_count = 0
            else:
                neg_term_count += 1
        if pos_term_count < 3:
            impure_sum_new, parafac_factors_new, boundaries_new =\
                offset_opt_func(offset_pos, impure_peak, quali_comp_db,
                                show_parafac_analytics)
            if impure_sum_new > impure_sum_opt:
                impure_sum_opt = impure_sum_new
                offset_opt = offset_pos
                parafac_factors_cur = parafac_factors_new
                boundaries_cur = boundaries_new
                pos_term_count = 0
            else:
                pos_term_count += 1
        
        # termination if three points in one direction without getting better
        if neg_term_count >= 3 and pos_term_count >= 3:
            termination_criterion = True
        
        # termination if algorithm runs to the constraints. Init value returned
        if offset > relative_distance_thresh * len(impure_peak.dataset.time):
            termination_criterion = True
            parafac_factors_cur, boundaries_cur, comp_tensor_shape, y_offset =\
                parafac(impure_peak, quali_comp_db, offset_init, show_parafac_analytics)

            offset_opt = offset_init

    return parafac_factors_cur, boundaries_cur, offset_opt, y_offset

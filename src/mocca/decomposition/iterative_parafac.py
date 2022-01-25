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


def get_all_comp_sum(parafac_factors, comp_tensor_shape,
                         show_parafac_analytics):
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


def iterative_parafac(impure_peak, quali_comp_db, show_parafac_analytics):
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
    iter_offset_cur = -impure_peak.offset
    # comp_tensor_shape should not change over iterative parafac runs
    parafac_factors_cur, boundaries_cur, comp_tensor_shape, y_offset =\
        parafac(impure_peak, quali_comp_db, iter_offset_cur, show_parafac_analytics)
    non_comp_sum_cur = get_all_non_comp_sum(parafac_factors_cur,
                                            comp_tensor_shape,
                                            show_parafac_analytics)
    comp_sum_cur = get_all_comp_sum(parafac_factors_cur,
                                    comp_tensor_shape,
                                    show_parafac_analytics)
    
    # first move positive
    iter_offset_new = iter_offset_cur + 1
    parafac_factors_new, boundaries_new, *_ = parafac(impure_peak,
                                                      quali_comp_db,
                                                      iter_offset_new,
                                                      show_parafac_analytics)
    non_comp_sum_new = get_all_non_comp_sum(parafac_factors_new,
                                            comp_tensor_shape,
                                            show_parafac_analytics)
    comp_sum_new = get_all_comp_sum(parafac_factors_new,
                                    comp_tensor_shape,
                                    show_parafac_analytics)

    # check if it gets better in positive direction
    if non_comp_sum_new < non_comp_sum_cur and comp_sum_new > comp_sum_cur:
        iter_max = int(math.ceil(0.05 * parafac_factors_new[1].shape[0])) + iter_offset_cur
        # find local minimum
        while (non_comp_sum_new < non_comp_sum_cur and
               comp_sum_new > comp_sum_cur and iter_offset_new < iter_max):
            non_comp_sum_cur = non_comp_sum_new
            comp_sum_cur = comp_sum_new
            parafac_factors_cur = parafac_factors_new
            boundaries_cur = boundaries_new
            iter_offset_cur = iter_offset_new
            iter_offset_new += 1
            parafac_factors_new, boundaries_new, *_ = parafac(impure_peak,
                                                              quali_comp_db,
                                                              iter_offset_new,
                                                              show_parafac_analytics)
            non_comp_sum_new = get_all_non_comp_sum(parafac_factors_new,
                                                    comp_tensor_shape,
                                                    show_parafac_analytics)
            comp_sum_new = get_all_comp_sum(parafac_factors_new,
                                            comp_tensor_shape,
                                            show_parafac_analytics)
    
    # if positive was wrong direction, let's go negative
    else:
        iter_offset_new = iter_offset_cur - 1
        parafac_factors_new, boundaries_new, *_ = parafac(impure_peak,
                                                         quali_comp_db,
                                                         iter_offset_new,
                                                         show_parafac_analytics)
        non_comp_sum_new = get_all_non_comp_sum(parafac_factors_new,
                                                comp_tensor_shape,
                                                show_parafac_analytics)
        comp_sum_new = get_all_comp_sum(parafac_factors_new,
                                        comp_tensor_shape,
                                        show_parafac_analytics)
        
        # check if it gets better in negative direction
        if non_comp_sum_new < non_comp_sum_cur and comp_sum_new > comp_sum_cur:
            iter_min = -int(math.floor(0.05 * parafac_factors_new[1].shape[0])) + iter_offset_cur
            # find local minimum
            while (non_comp_sum_new < non_comp_sum_cur and
                   comp_sum_new > comp_sum_cur and iter_offset_new > iter_min):
                non_comp_sum_cur = non_comp_sum_new
                comp_sum_cur = comp_sum_new
                parafac_factors_cur = parafac_factors_new
                boundaries_cur = boundaries_new
                iter_offset_cur = iter_offset_new
                iter_offset_new -= 1
                
                parafac_factors_new, boundaries_new, *_ = parafac(impure_peak,
                                                                 quali_comp_db,
                                                                 iter_offset_new,
                                                                 show_parafac_analytics)
                non_comp_sum_new = get_all_non_comp_sum(parafac_factors_new,
                                                        comp_tensor_shape,
                                                        show_parafac_analytics)
                comp_sum_new = get_all_comp_sum(parafac_factors_new,
                                                comp_tensor_shape,
                                                show_parafac_analytics)
        
        # if negative also does not give progress, start point is minimum
    return parafac_factors_cur, boundaries_cur, iter_offset_cur, y_offset

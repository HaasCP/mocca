#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 09:17:34 2022

@author: haascp
"""


def get_non_comp_sum(parafac_factors, start_slice, end_slice):
    """
    Get summed integrals from all components except the highest inegral component
    in one slice of the PARAFAC factors.
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
    Get summed integrals from the highest inegral component in one slice of
    PARAFAC factors.
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
    """
    Iterative PARAFAC approach is a maximization problem. Based on the component
    tensor data shape, the maximal sum of integrals is found for every component
    part in the data tensor. The integrals of these component integrals are
    summed up. This sum is the maximization objective for the iterative PARAFAC
    approach.
    """
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
    """
    Iterative PARAFAC approach is a maximization problem. This objective function
    assumes that the PARAFAC model is best, when the overall integral over all
    slices is maximized.
    """
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

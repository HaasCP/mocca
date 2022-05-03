# -*- coding: utf-8 -*-
"""
Created on Mon May  2 13:45:47 2022

@author: CPH
"""
import numpy as np

from mocca.dad_data.utils import trim_data
from mocca.dad_data.process_gradientdata import bsl_als

def get_bsl_corrected_comp_data(comp_data):
    """
    Returns baseline-corrected compound data which will be used to be compared
    against the baseline-corrected blank gradient data.
    """
    comp_data_bsl = bsl_als(comp_data)
    return comp_data_bsl


def get_mse_comp_grad(comp_data, grad_data):
    """
    Returns mean squared error between two arrays.
    """
    difference_array = np.subtract(comp_data, grad_data)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    return mse


def get_thresh_idxs(dataset, peaks_high_pass, peaks_low_pass):
    times = dataset.time
    if peaks_high_pass:
        high_pass_idx = next((times.index(n) for n in times if
                              n > peaks_high_pass), len(times))
    else:
        high_pass_idx = 0

    if peaks_low_pass:
        low_pass_idx = next((times.index(n) for n in times if
                              n > peaks_low_pass), len(times))
    else:
        low_pass_idx = len(times)

    return high_pass_idx, low_pass_idx


def get_opt_grad_data(comp_dataset, grad_dataset, peaks_high_pass, peaks_low_pass,
                      relative_distance_thresh):
    """
    
    """
    comp_data_bsl = get_bsl_corrected_comp_data(comp_dataset.data)

    high_pass_idx, low_pass_idx = get_thresh_idxs(comp_dataset, peaks_high_pass,
                                                  peaks_low_pass)
    grad_data = grad_dataset.data[:, high_pass_idx:low_pass_idx]
    comp_data_bsl = comp_data_bsl[:, high_pass_idx:low_pass_idx]

    len_iterator = int(relative_distance_thresh / 2 * comp_data_bsl.shape[1])
    offset_iterator = [i - len_iterator for i in
                       list(range(len_iterator * 2 + 1))]

    mse_opt = float("inf")

    iter_objective_func = []
    for offset in offset_iterator:
        pass

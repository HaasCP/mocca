# -*- coding: utf-8 -*-
"""
Created on Mon May  2 13:45:47 2022

@author: CPH
"""
import numpy as np

from mocca.dad_data.utils import trim_data
from mocca.dad_data.process_gradientdata import bsl_als
import matplotlib.pyplot as plt

def get_bsl_corrected_comp_data(comp_data):
    """
    Returns baseline-corrected compound data which will be used to be compared
    against the baseline-corrected blank gradient data.
    """
    comp_data_bsl = bsl_als(comp_data)
    return comp_data_bsl


def get_thresh_idxs(dataset, peaks_high_pass, peaks_low_pass):
    """
    Get indices in the time domain which represent the high pass and low pass
    filters.
    """
    times = list(dataset.time)
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


def get_mse_comp_grad(comp_data, grad_data, high_pass_idx, low_pass_idx):
    """
    Returns item-wise mean squared error of two arrays in between the given
    filter indices.
    """
    min_length = min(comp_data.shape[1], grad_data.shape[1], low_pass_idx)
    comp_data = comp_data[:, high_pass_idx:min_length]
    grad_data = grad_data[:, high_pass_idx:min_length]

    difference_array = np.subtract(comp_data, grad_data)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    return mse


def get_iter_grad_data(grad_data, i, comp_data):
    """
    Returns a shifted gradient data array.
    """
    if i < 0:
        grad_data = grad_data[:, -i:]  # left shiftes
    elif i > 0:
        grad_data = np.c_[comp_data[:, :i], grad_data]
    return grad_data




def get_opt_grad_data(comp_dataset, grad_dataset, peaks_high_pass, peaks_low_pass,
                      relative_distance_thresh):
    """
    
    """
    comp_data_bsl = get_bsl_corrected_comp_data(comp_dataset.data)
    grad_data = grad_dataset.data

    high_pass_idx, low_pass_idx = get_thresh_idxs(comp_dataset, peaks_high_pass,
                                                  peaks_low_pass)
    
    comp_data_bsl = comp_data_bsl[:, high_pass_idx:low_pass_idx]

    len_iterator = int(relative_distance_thresh / 2 * comp_data_bsl.shape[1])
    offset_iterator = [i - len_iterator for i in
                       list(range(len_iterator * 2 + 1))]

    mse_opt = float("inf")
    grad_data_opt = None

    iter_objective_func = []
    for offset in offset_iterator:
        grad_data_iter = get_iter_grad_data(grad_data, offset, comp_data_bsl)
        mse = get_mse_comp_grad(comp_data_bsl, grad_data_iter, high_pass_idx,
                                low_pass_idx)
        print(offset)
        print(mse)
        print(grad_data_iter[:3, :10])
        plt.contour(grad_data_iter[:, :50])
        plt.show()
        if mse < mse_opt:
            mse_opt = mse
            grad_data_opt = grad_data_iter
    return grad_data_opt
        



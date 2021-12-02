#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 12:16:35 2021

@author: haascp
"""
import numpy as np
import matplotlib.pyplot as plt

from mocca.peak.process.purity_funcs import (get_trimmed_peak_data, get_max_loc,
                                     get_noise_variance, get_correls,
                                     get_agilent_thresholds,
                                     get_purity_value_agilent,
                                     predict_purity_unimodal,
                                     get_pca_explained_variance)

def predict_peak_purity(peak, param=2.5, show_analytics=False):
    """
    Returns peak purity prediction by performing the described test sequence.
    Plots and prints infromation about the peak purity prediction.
    """
    peak_data = get_trimmed_peak_data(peak)
    max_loc = get_max_loc(peak_data)
    noise_variance = get_noise_variance(peak)
    correls = get_correls(peak_data, max_loc)    
    agilent_thresholds = get_agilent_thresholds(peak_data, max_loc, 
                                                noise_variance, param)
    
    test_agilent = get_purity_value_agilent(peak_data, correls, 
                                            agilent_thresholds)
    test_unimodality = predict_purity_unimodal(correls)
    test_pca = get_pca_explained_variance(peak_data)
    test_correls_1 = np.min(correls)
    test_correls_2 = np.mean(correls)
    
    
    # peak purity algorithm
    
    #  if agilent threshold reached, then probably pure
    if test_agilent > 0.9:
        return True
    #  for pure peak, correlation array emperically expected to be unimodal
    if not test_unimodality:
        return False
    #  if pca big enough, then probably pure
    if test_pca > 0.995:
        return True
    #  if any correlation is < 0.9, then probably impure somewhere
    if test_correls_1 < 0.9:
        return False
    #  otherwise, check if correlation shows that it is reasonably pure
    if test_correls_1 > 0.95:
        return True
    # if peak is pure, average correlation across all peaks should be high
    if test_correls_2 > 0.98:
        return True
    return False

    if show_analytics:
        plt.plot(correls[0])
        plt.plot(agilent_thresholds)
        plt.show()
        plt.plot(correls[0])
        plt.plot(correls[1])
        plt.show()
        for i in range(peak_data.shape[1]):
            plt.plot(peak_data[:, i])
        plt.show()
        print(f"Agilent Threshold (True for >0.9): {test_agilent} \n"
              f"Unimodality Test (False for False): {test_unimodality} \n"
              f"PCA Variance Explained (True for >0.995): {test_pca} \n"
              f"Minimum Correlation (False for <0.9): {test_correls_1} \n"
              f"Minimum Correlation (True for >0.95): {test_correls_1} \n"
              f"Average Correlation (True for >0.98): {test_correls_2} \n")

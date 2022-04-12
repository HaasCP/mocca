#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:55:25 2021

@author: haascp
"""

from mocca.peak.utils import get_peak_data, is_unimodal

import numpy as np
from sklearn.decomposition import PCA


def get_trimmed_peak_data(peak):
    """
    Returns peak data trimmed with cut edges of the peak to 5% of max absorbance
    to avoid noise artifacts.
    """
    peak_data = get_peak_data(peak)
    return peak_data[:, np.sum(peak_data, axis=0) > 0.05 *
                     np.max(np.sum(peak_data, axis=0))]


def get_max_loc(peak_data):
    """
    Returns the maximum location of the given peak data.
    """
    return np.argmax(np.sum(peak_data, axis=0))


def get_noise_variance(peak):
    """
    Filters dataset with only timepoints whose max absorbance at
    any wavelength is below 1% of max absorbance. Returns the average of the
    variance over all wavelengths.
    """
    noise_data = peak.dataset.data[:, np.max(peak.dataset.data, axis=0) <
                                   0.01 * np.max(peak.dataset.data)]
    return np.mean(np.var(noise_data, axis=0))


def get_correls(peak_data, max_loc):
    """
    Get a list with correlation coefficients of UV-Vis spectra at every
    timepoint with reference to the UV-Vis spectrum at maximum absorbance.
    """
    correls_to_max = [(np.corrcoef(peak_data[:, i],
                                   peak_data[:, max_loc])[0, 1])**2
                      for i in range(peak_data.shape[1])]
    return correls_to_max


def get_agilent_thresholds(peak_data, max_loc, noise_variance, param=2.5):
    """
    Returns the thresholds calculated by the Agilent purity algorithm.
    """
    agilent_thresholds = [(max(0, 1 - param *
                               (noise_variance / np.var(peak_data[:, i]) +
                                noise_variance / np.var(peak_data[:, max_loc]))))**2  # noqa: E501
                          for i in range(peak_data.shape[1])]
    return agilent_thresholds


def get_purity_value_agilent(peak_data, correls, agilent_thresholds):
    """
    Uses Agilent's peak purity algorithm to predict purity of peak. Param
    gives strictness of test (original was 0.5, which is more strict)
    """
    # check if > 90% of the points are greater than the modified agilent threshold.
    agilent_test = np.sum(np.greater(correls,
                                     agilent_thresholds)) / peak_data.shape[1]
    return agilent_test


def predict_purity_unimodal(correls):
    """
    Checks for unimodality of a peak by an averaging filter of length 3
    on the correlation vector to the maximum
    https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-numpy-scipy
    """
    return is_unimodal(np.convolve(correls, np.ones(3), 'valid') / 3, 0.999)


def get_pca_explained_variance(peak_data):
    """
    Calculates the ration of explained variance by the first principal
    component of the devonvoluted peak data.
    """
    pca = PCA(n_components=1)
    pca.fit(peak_data)
    return pca.explained_variance_ratio_[0]

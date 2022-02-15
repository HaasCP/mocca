#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 08:57:14 2022

@author: haascp
"""
from sklearn.decomposition import PCA
from tensorly.decomposition import non_negative_parafac_hals
import matplotlib.pyplot as plt 

from mocca.peak.utils import get_peak_data
from mocca.decomposition.utils import normalize_parafac_factors
from mocca.decomposition.data_tensor import get_parafac_tensor


def estimate_num_components_pca(data_tensor, impure_peak, show_parafac_analytics):
    """
    Returns an estimate of the number of components in the impure peak.
    """
    pca_data = data_tensor.reshape(get_peak_data(impure_peak).shape[0], -1)
    
    pca = PCA(n_components=10)
    _ = pca.fit_transform(pca_data)

    threshold = 0.9995
    cum_sum = 0
    for idx in range(10):
        cum_sum += pca.explained_variance_ratio_[idx]
        if show_parafac_analytics:
            print(f"Cumulative variance after {idx + 1} PCA components: {cum_sum}")
        if cum_sum > threshold:
            return idx + 1


def print_parafac_analytics(normalized_spectra, normalized_elution, normalized_concs,
                      boundaries, relevant_comp):
    """
    Prints out PARAFAC decomposition model results. Used for debugging and dev.
    """
    print(f"PARAFAC data tensor has boundaries in the time axis of {boundaries}.")
    plt.plot(normalized_spectra)
    print(f"Compound used for PARAFAC data tensor: {relevant_comp.compound_id}")
    plt.plot([val / (max(relevant_comp.spectrum) / normalized_spectra.max()) 
              for val in relevant_comp.spectrum], "--")
    plt.show()
    plt.plot(normalized_elution)
    plt.show()
    plt.plot(normalized_concs)
    plt.show()


def parafac(impure_peak, quali_comp_db, iter_offset, show_parafac_analytics):
    """
    PARAFAC decomposition processing routine of impure peaks. An iteration offset
    is introduced to allow for iterative PARAFAC approach.
    """
    if show_parafac_analytics:
        print(f"----- new PARAFAC decomposition with iteration offset {iter_offset}"
              " -----")
    data_tensor, boundaries, relevant_comp, comp_tensor_shape, y_offset =\
        get_parafac_tensor(impure_peak, quali_comp_db, iter_offset,
                           show_parafac_analytics)

    pca_n_comps = estimate_num_components_pca(data_tensor, impure_peak,
                                              show_parafac_analytics)
    n_comps = max(pca_n_comps, 2)
    
    weights, factors = non_negative_parafac_hals(data_tensor, rank=n_comps,
                                                 init='svd', n_iter_max=1000,
                                                 verbose = False, tol=1e-9)

    spectra = factors[0]
    elutions = factors[1]
    integrals = factors[2]

    normalized_spectra, normalized_elution, normalized_integrals =\
        normalize_parafac_factors(spectra, elutions, integrals)

    if show_parafac_analytics:
        print(f"Estimated number of components is {n_comps}")
        print_parafac_analytics(normalized_spectra, normalized_elution,
                                normalized_integrals, boundaries, relevant_comp)
        print(f"integral_array = {normalized_integrals}")

    parafac_factors = (normalized_spectra, normalized_elution, normalized_integrals)
    
    return parafac_factors, boundaries, comp_tensor_shape, y_offset

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


def estimate_num_components_pca(data_tensor, impure_peak):
    """
    Returns an estimate of the number of components in the impure peak.
    """
    pca_data = data_tensor.reshape(get_peak_data(impure_peak).shape[0], -1)
    
    pca = PCA(n_components=10)
    _ = pca.fit_transform(pca_data)

    threshold = 0.998
    cum_sum = 0
    for idx in range(10):
        cum_sum += pca.explained_variance_ratio_[idx]
        if cum_sum > threshold:
            return idx + 1


def print_parafac_analytics(normalized_spectra, normalized_elution, normalized_concs,
                      boundaries, relevant_comps):
    """
    Prints out PARAFAC decomposition model results. Used for debugging and dev.
    """
    print(f"PARAFAC data tensor has boundaries in the time axis of {boundaries}.")
    plt.plot(normalized_spectra)
    for comp in relevant_comps:
        print(f"Compound used for PARAFAC data tensor: {comp.compound_id}")
        plt.plot([val / (max(comp.spectrum) / normalized_spectra.max()) 
                  for val in comp.spectrum], "--")
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
    print(f"----- new PARAFAC decomposition with iteration offset {iter_offset}"
          " -----")
    data_tensor, boundaries, relevant_comps, comp_tensor_shape =\
        get_parafac_tensor(impure_peak, quali_comp_db, iter_offset)

    pca_n_comps = estimate_num_components_pca(data_tensor, impure_peak)
    n_comps = max(len(relevant_comps) + 1, pca_n_comps)
    
    weights, factors = non_negative_parafac_hals(data_tensor, rank=n_comps,
                                                 init='svd', n_iter_max=1000,
                                                 verbose = False, tol=1e-9)

    spectra = factors[0]
    elutions = factors[1]
    integrals = factors[2]

    normalized_spectra, normalized_elution, normalized_integrals =\
        normalize_parafac_factors(spectra, elutions, integrals)

    if show_parafac_analytics:
        print_parafac_analytics(normalized_spectra, normalized_elution,
                                normalized_integrals, boundaries, relevant_comps)
        print(f"integral_array = {normalized_integrals}")

    parafac_factors = (normalized_spectra, normalized_elution, normalized_integrals)
    
    return parafac_factors, boundaries, comp_tensor_shape


#print(len(factors[1]))
#for data in factors[1]:
#    print(data.shape)
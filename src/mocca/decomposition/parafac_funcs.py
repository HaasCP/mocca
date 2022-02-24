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
from mocca.decomposition.data_tensor import get_parafac_tensor
from mocca.decomposition.model import ParafacModel


def estimate_pca_n_comps(data_tensor, impure_peak, show_parafac_analytics):
    """
    Returns an estimate of the number of components in the impure peak.
    """
    pca_data = data_tensor.reshape(get_peak_data(impure_peak).shape[0], -1)

    pca = PCA(n_components=10)
    _ = pca.fit_transform(pca_data)

    #  same thresh than in peak.check, check_peak_purity
    threshold = 0.9995
    cum_sum = 0
    pca_explained_variance = []
    for idx in range(10):
        cum_sum += pca.explained_variance_ratio_[idx]
        if show_parafac_analytics:
            print(f"Cumulative variance after {idx + 1} PCA components: {cum_sum}")
        pca_explained_variance.append(cum_sum)
        if cum_sum > threshold:
            return idx + 1, pca_explained_variance


def print_parafac_analytics(parafac_model):
    """
    Prints out PARAFAC decomposition model results. Used for debugging and dev.
    """
    tensor = parafac_model.data_tensor
    print("PARAFAC data tensor has boundaries in the time axis of "
          f"{tensor.boundaries}.")
    normalized_spectra, normalized_elution, normalized_integrals =\
        parafac_model.factors
    print(f"Estimated number of components is {parafac_model.n_comps}")
    plt.plot(normalized_spectra)
    print("Compound used for PARAFAC data tensor: "
          f"{tensor.relevant_comp.compound_id}")
    plt.plot([val / (max(tensor.relevant_comp.spectrum) / normalized_spectra.max())
              for val in tensor.relevant_comp.spectrum], "--")
    plt.show()
    plt.plot(normalized_elution)
    plt.show()
    plt.plot(normalized_integrals)
    plt.show()
    print(f"integral_array = {normalized_integrals}")


def parafac(impure_peak, quali_comp_db, iter_offset, show_parafac_analytics):
    """
    PARAFAC decomposition processing routine of impure peaks. An iteration offset
    is introduced to allow for iterative PARAFAC approach.
    """
    if show_parafac_analytics:
        print(f"----- new PARAFAC decomposition with iteration offset {iter_offset}"
              " -----")
    data_tensor = get_parafac_tensor(impure_peak, quali_comp_db, iter_offset,
                                     show_parafac_analytics)

    pca_n_comps, pca_explained_variance = estimate_pca_n_comps(data_tensor.tensor,
                                                               impure_peak,
                                                               show_parafac_analytics)
    n_comps = max(pca_n_comps, 2)

    weights, factors = non_negative_parafac_hals(data_tensor.tensor, rank=n_comps,
                                                 init='svd', n_iter_max=1000,
                                                 verbose=False, tol=1e-9)

    parafac_model = ParafacModel(impure_peak, n_comps, pca_explained_variance,
                                 weights, factors, data_tensor, iter_offset)

    if show_parafac_analytics:
        print_parafac_analytics(parafac_model)

    return parafac_model

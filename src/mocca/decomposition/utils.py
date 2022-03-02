#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 09:01:51 2022

@author: haascp
"""
import numpy as np
import itertools

from mocca.dad_data.utils import sum_absorbance_by_time


def check_comp_overlap(peak, comp):
    """
    Checks if a given peak overlaps with a given component.
    """
    return comp.left <= peak.left - peak.offset <= comp.right\
        or peak.left - peak.offset <= comp.left <= peak.right - peak.offset


def check_any_compound_overlap(peak, quali_comp_db):
    """
    Checks if a given peak overlaps with any component in the quali_comp_db.
    """
    return any((check_comp_overlap(peak, comp) and 'unknown' not in comp.compound_id
                and 'impurity' not in comp.compound_id) for comp in quali_comp_db)


def check_same_uvvis(parafac_model, spectrum_correl_coef_thresh):
    """
    Checks if any two parafac components share the same UV-Vis trace.
    """
    spectra = parafac_model.factors[0]
    if all(np.corrcoef(spectra[:, i], spectra[:, j])[0, 1] > 
           spectrum_correl_coef_thresh for i, j in 
           itertools.combinations(list(range(parafac_model.n_comps)), 2)):
        return True
    else:
        return False


def check_comp_in_impure(parafac_model, absorbance_threshold):
    """
    Checks if the maximum absorbance of the known PARAFAC component in the impure
    peak exceeds the absorbance threshold.
    """
    integrals = parafac_model.factors[2]
    comp_idx = []
    comp_integrals = []
    for integral_slice in integrals[:-1]:
        comp_integrals.append(integral_slice.max())
        comp_idx.append(integral_slice.argmax())
    if not all(idx == comp_idx[0] for idx in comp_idx):
        return True
    else:
        comp_idx = comp_idx[0]
    spectrum_max = parafac_model.factors[0][-1][comp_idx].max()
    elution_max = parafac_model.factors[1][-1][comp_idx].max()
    integral = integrals[-1][comp_idx].max()
    abs_max_comp_in_impure = spectrum_max * elution_max * integral
    if abs_max_comp_in_impure < absorbance_threshold:
        return False
    else: 
        return True


def check_absorbance_thresh(parafac_peak, absorbance_threshold):
    """
    Checks if maximum absorbance in synthetically created PARAFAC peak dataset
    exceeds absorbance threshold.
    """
    max_absorbance = np.max(sum_absorbance_by_time(parafac_peak.dataset.data))
    return max_absorbance > absorbance_threshold

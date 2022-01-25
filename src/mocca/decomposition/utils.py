#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 09:01:51 2022

@author: haascp
"""
import numpy as np

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
    return any((check_comp_overlap(peak, comp) and not 'unknown' in comp.compound_id
                and not 'impurity' in comp.compound_id) for comp in quali_comp_db)


def normalize_parafac_factors(spectra, elutions, integrals):
    """
    Normalizes the calculated PARAFAC factors to bring them in the format of
    the original data.
    """
    spectral_norm_val = np.sum(spectra, axis=0)
    elution_norm_val = np.sum(elutions, axis=0)

    normalized_spectra = spectra / spectral_norm_val
    normalized_elution = elutions / elution_norm_val
    normalized_integrals = integrals * spectral_norm_val * elution_norm_val

    return normalized_spectra, normalized_elution, normalized_integrals

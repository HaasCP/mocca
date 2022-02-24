#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 11:15:42 2022

@author: haascp
"""
import numpy as np
from dataclasses import dataclass
from typing import Union

from mocca.components.models import QualiComponent
from mocca.peak.models import CorrectedPeak, IntegratedPeak


@dataclass(frozen=True)
class DataTensor():
    """
    Model of data tensors used as input for the PARAFAC decomposition algorithm.
    """
    tensor : np.ndarray
    boundaries : tuple
    relevant_comp : QualiComponent
    comp_tensor_shape : tuple
    y_offset : float


@dataclass()
class ParafacModel():
    """
    Stores all relevant information of a PARAFAC model.
    """
    impure_peak : Union[CorrectedPeak, IntegratedPeak]
    n_comps : int
    pca_explained_variance : list
    weights : list
    factors : list
    data_tensor : DataTensor

    def __post_init__(self):
        self._normalize_factors()

    def _normalize_factors(self):
        """
        Normalizes the calculated PARAFAC factors to bring them in the format of
        the original data.
        """
        spectra = self.factors[0]
        elutions = self.factors[1]
        integrals = self.factors[2]

        spectral_norm_val = np.sum(spectra, axis=0)
        elution_norm_val = np.sum(elutions, axis=0)

        normalized_spectra = spectra / spectral_norm_val
        normalized_elution = elutions / elution_norm_val
        normalized_integrals = integrals * spectral_norm_val * elution_norm_val

        self.factors = [normalized_spectra, normalized_elution, normalized_integrals]

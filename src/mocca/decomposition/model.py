#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 11:15:42 2022

@author: haascp
"""
import numpy as np
from dataclasses import dataclass, field
from typing import Union, List, Optional

from mocca.components.models import QualiComponent
from mocca.peak.models import CorrectedPeak, IntegratedPeak

from mocca.decomposition.utils import (check_same_uvvis, check_summed_factor_uvvis,
                                       check_comp_in_impure)
from mocca.peak.resolve_impure import create_pure_peak, create_parafac_peak


@dataclass()
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
    iter_offset : int
    iter_objective_func : list = None
    peaks : Optional[List[Union[CorrectedPeak, IntegratedPeak]]] = None
    impure_mse : float = field(init=False)

    def __post_init__(self):
        self._normalize_factors()
        self._calculate_impure_mse()

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

        normalized_spectra = np.divide(spectra, spectral_norm_val,
                                       out=np.zeros_like(spectra),
                                       where=spectral_norm_val != 0)
        normalized_elution = np.divide(elutions, elution_norm_val,
                                       out=np.zeros_like(elutions),
                                       where=elution_norm_val != 0)
        normalized_integrals = integrals * spectral_norm_val * elution_norm_val

        self.factors = [normalized_spectra, normalized_elution, normalized_integrals]

    def _calculate_impure_mse(self):
        """
        Calculates the mean square error of the summed PARAFAC components for the
        impure peak slice and the actual impure peak.
        """
        impure_data = self.data_tensor.tensor[:, :, -1]
        parafac_data = np.zeros_like(impure_data)
        for comp_i in range(self.n_comps):
            n_wls = len(self.factors[0][:, comp_i])
            comp_wl = self.factors[0][:, comp_i].reshape(n_wls, 1)
            comp_time = self.factors[1][:, comp_i]
            comp_integral = self.factors[2][:, comp_i][-1]
            comp_data = comp_wl * comp_time * comp_integral
            parafac_data = np.add(parafac_data, comp_data)
        difference_array = np.subtract(impure_data, parafac_data)
        squared_array = np.square(difference_array)
        mse = squared_array.mean()
        self.impure_mse = mse

    def create_parafac_peaks(self, absorbance_threshold,
                             spectrum_correl_coef_thresh):
        """
        If two UV-Vis traces in the PARAFAC components are too similar, no
        PARAFAC peaks are created. PARAFAC peaks' dataseta are created synthetically
        generated by using the PARAFAC factors of the model and filling the rest
        of the array up with zeros. PARAFAC peaks get an index of -impure_peak.idx
        """
        if (check_same_uvvis(self, spectrum_correl_coef_thresh) or
                check_summed_factor_uvvis(self, spectrum_correl_coef_thresh)):
            parafac_peaks = [create_pure_peak(self.impure_peak)]
        elif not check_comp_in_impure(self, absorbance_threshold):
            parafac_peaks = []
        else:
            parafac_peaks = []
            for i in range(self.n_comps):
                #  get factors for one parafac comonent
                parafac_peak = create_parafac_peak(i, self)
                parafac_peaks.append(parafac_peak)
        self.peaks = parafac_peaks

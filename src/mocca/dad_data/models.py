#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 13:16:51 2021

@author: haascp
"""

from dataclasses import dataclass, field, InitVar
from typing import List
import numpy as np
import copy

from mocca.dad_data.utils import trim_data
from mocca.dad_data.process_gradientdata import bsl_als
from mocca.dad_data.apis.chemstation import read_chemstation
from mocca.dad_data.apis.labsolutions import read_labsolutions
from mocca.dad_data.apis.custom import read_custom_data

import mocca.peak.models


@dataclass()
class DadData():
    """
    Base class for HPLC-DAD data.
    """
    hplc_system_tag : str
    experiment : InitVar["mocca.user_interaction.user_objects.HplcInput"]
    wl_high_pass : InitVar[float] = None
    wl_low_pass : InitVar[float] = None

    # set during initialization
    path : str = field(init=False)
    data : np.ndarray = field(init=False)
    time : np.ndarray = field(init=False)
    wavelength : np.ndarray = field(init=False)
    warnings : List[str] = field(init=False)

    def __post_init__(self, experiment, wl_high_pass, wl_low_pass):
        """
        Process init variables to set attributes.
        """
        self.warnings = []
        self._set_path(experiment)
        self._read_data(wl_high_pass, wl_low_pass, experiment)
        experiment.processed = True

    def __eq__(self, other):
        """
        Checks if two DAD data containers are equal by comparing its path and data.
        """
        if not isinstance(other, type(self)):
            # don't attempt to compare against unrelated types
            return False
        return (self.path == other.path and
                (self.data.__array_interface__['data'] ==
                 other.data.__array_interface__['data']))

    def _set_path(self, experiment):
        """
        Sets path to the data file.
        """
        self.path = experiment.path

    def _read_data(self, wl_high_pass, wl_low_pass, experiment):
        """
        Read the data file, preprocess it and filter it in order to obtain
        standardized data format for all HPLC systems.
        """
        if self.hplc_system_tag == 'chemstation':
            data, time, wavelength = read_chemstation(self.path, wl_high_pass,
                                                      wl_low_pass)
        elif self.hplc_system_tag == 'labsolutions':
            data, time, wavelength = read_labsolutions(self.path, wl_high_pass,
                                                       wl_low_pass)
        elif self.hplc_system_tag == 'custom':
            data, time, wavelength = read_custom_data(experiment)
        else:
            raise ValueError(f"Given hplc_system_tag  {self.hplc_system_tag} is "
                             "not defined in MOCCA.")
        self.data = data
        self.time = time
        self.wavelength = wavelength


@dataclass(eq=False)
class GradientData(DadData):
    """
    Data container for gradient HPLC-DAD data.
    """
    original_data : np.ndarray = field(init=False)

    def __post_init__(self, experiment, wl_high_pass, wl_low_pass):
        """
        Process and baseline-correct given data.
        """
        super().__post_init__(experiment, wl_high_pass, wl_low_pass)
        self.original_data = copy.deepcopy(self.data)
        self.data = bsl_als(self.data)


@dataclass(eq=False)
class CompoundData(DadData):
    """
    Data container for HPLC-DAD data with peaks originating from compounds.
    """
    def __post_init__(self, experiment, wl_high_pass, wl_low_pass):
        """
        Baseline-corrects the given HPLC-DAD data.
        """
        super().__post_init__(experiment, wl_high_pass, wl_low_pass)
        if experiment.gradient is not None:
            self._subtract_baseline(experiment.gradient.dataset)

    def _trim_data(self, length):
        """Trims the data in the time dimension to the length provided"""
        self.data, self.time = trim_data(data=self.data, time=self.time, length=length)

    def _subtract_baseline(self, gradient):
        """Subtracts the baseline of the gradient numpy array from self.data."""
        self._trim_data(gradient.data.shape[1])
        self.data = self.data - gradient.data[:, :self.data.shape[1]]


@dataclass
class ParafacData():
    """
    Data container for synthetic data generated from PARAFAC models.
    """
    # https://www.python.org/dev/peps/pep-0484/#forward-references
    impure_peak : InitVar['mocca.peak.models.CorrectedPeak']
    parafac_comp_tensor : InitVar[tuple]
    boundaries : InitVar[tuple]
    shift : InitVar[int]
    y_offset : InitVar[float]

    def __post_init__(self, impure_peak, parafac_comp_tensor, boundaries,
                      shift, y_offset):
        """
        Process init variables to set attributes and create a synthetic HPLC-DAD
        dataset.
        """
        # https://github.com/python/mypy/issues/9254
        self.hplc_system_tag = impure_peak.dataset.hplc_system_tag
        self.path = impure_peak.dataset.path
        self.time = impure_peak.dataset.time
        self.wavelength = impure_peak.dataset.wavelength
        self.data = np.zeros((len(self.wavelength), len(self.time)))
        self._make_data_from_parafac_peak(impure_peak, parafac_comp_tensor,
                                          boundaries, shift, y_offset)

    def _make_data_from_parafac_peak(self, impure_peak, parafac_comp_tensor,
                                     boundaries, shift, y_offset):
        """
        Create a synthetic dataset generated from an impure peak and a
        corresponding PARAFAC model.
        """
        # make 2D data corresponding to parafac-generated spectra and elution
        spectrum = parafac_comp_tensor[0].reshape(len(self.wavelength), 1)
        retention = parafac_comp_tensor[1].reshape(1, boundaries[1] - boundaries[0] + 1)
        integral = parafac_comp_tensor[2][-1]  # reaction run is last in run dimension
        parafac_peak_data = spectrum * retention * integral

        # replace self data with parafac peak data
        left = boundaries[0] - shift  # impure_peak.offset already included in
        right = boundaries[1] - shift + 1
        self.data[:, left:right] = parafac_peak_data + y_offset
        self.data[self.data < 0] = 0

    def __eq__(self, other):
        """
        Checks if two ParafacData objects are the same based on the path and the
        data.
        """
        if not isinstance(other, type(self)):
            # don't attempt to compare against unrelated types
            return False
        return (self.path == other.path and
                (self.data.__array_interface__['data'] ==
                 other.data.__array_interface__['data']))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 13:16:51 2021

@author: haascp
"""

from dataclasses import dataclass, field, InitVar
from typing import Dict, List, Optional

import numpy as np
from mocca.dad_data.utils import absorbance_to_array, apply_filter, trim_data
from mocca.dad_data.process_gradientdata import bsl_als

from mocca.dad_data.apis.chemstation_api import read_csv_agilent, tidy_df_agilent
from mocca.dad_data.apis.labsolutions_api import read_txt_shimadzu

from mocca.campaign.experiment import Experiment

import mocca.peak.models

# Parameter inheritance issues solved as shown in:
# https://stackoverflow.com/questions/51575931/class-inheritance-in-python-3-7-dataclasses
@dataclass()
class _DadDataBase():
    hplc_system_tag : str
    # input information, most importantly the path
    experiment : Experiment


@dataclass()
class _DadDataDefaultsBase():
    wl_high_pass : InitVar[float] = None
    wl_low_pass : InitVar[float] = None
    detector_limit : InitVar[float] = None


@dataclass()
class DadData(_DadDataDefaultsBase, _DadDataBase):
    # set during initialization
    path : str = field(init=False)
    detector_limit : int = field(init=False)
    data : np.ndarray = field(init=False)
    time : np.ndarray = field(init=False)
    wavelength : np.ndarray = field(init=False)

    def __post_init__(self, wl_high_pass, wl_low_pass, detector_limit):
        self._set_path()
        self._set_detector_limit(detector_limit)
        self._read_data(wl_high_pass, wl_low_pass)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            # don't attempt to compare against unrelated types
            raise ValueError("Both DAD datasets must be of the same type!")
        return self.path == other.path
    
    def _set_path(self):
        self.path = self.experiment.path

    def _set_detector_limit(self, detector_limit):
        if detector_limit is None:
            if self.hplc_system_tag == 'chemstation':
                self.detector_limit = 2000
            elif self.hplc_system_tag == 'labsolutions':
                self.detector_limit = 2000
            else:
                raise AttributeError("HPLC System Tag {} not supported!".format(self.hplc_system_tag))
        else:
            self.detector_limit = detector_limit

    def _read_data(self, wl_high_pass, wl_low_pass):
        if self.hplc_system_tag == 'chemstation':
            df = read_csv_agilent(self.path)
            df = apply_filter(tidy_df_agilent(df, wl_high_pass, wl_low_pass))
        elif self.hplc_system_tag == 'labsolutions':
            df = read_txt_shimadzu(self.path, wl_high_pass, wl_low_pass)
            df = apply_filter(df)

        self.data = absorbance_to_array(df)
        self.time = df.time.unique()
        self.wavelength = df.wavelength.unique()

@dataclass 
class _ParafacDataBase(_DadDataBase):
    parafac_peak : InitVar['mocca.peak.models.ParafacPeak']  # https://www.python.org/dev/peps/pep-0484/#forward-references
    original_dataset : InitVar[DadData]

@dataclass
class ParafacData(DadData, _ParafacDataBase):
    def __post_init__(self, parafac_peak, original_dataset, wl_high_pass, wl_low_pass):
        # https://github.com/python/mypy/issues/9254
        self.detector_limit = np.inf
        self.time = original_dataset.time
        self.wavelength = original_dataset.wavelength
        self.data = np.zeros((len(self.wavelength), len(self.time)))
        self._make_data_from_parafac_peak(parafac_peak)

    def _make_data_from_parafac_peak(self, parafac_peak):
        # make 2D data corresponding to parafac-generated spectra and elution
        parafac_peak_data = parafac_peak.spectra.reshape(len(self.wavelength), 1) \
                          * parafac_peak.elution.reshape(1, parafac_peak.right + 1 - parafac_peak.left) \
                          * parafac_peak.integral
        # replace self data with parafac peak data
        self.data[:, parafac_peak.left + parafac_peak.offset:parafac_peak.right + parafac_peak.offset + 1] = parafac_peak_data


@dataclass(eq=False)
class GradientData(DadData):
    def __post_init__(self, wl_high_pass, wl_low_pass):
        super().__post_init__(wl_high_pass, wl_low_pass)
        self.data = bsl_als(self.data)


@dataclass()
class _CompoundDataBase(_DadDataBase):
    # gradient for baseline correction
    gradient : InitVar[GradientData]


@dataclass()
class _CompoundDataDefaultsBase(_DadDataDefaultsBase):
    # data properties to be set by class, do not initialize
    warnings : List[str] = field(default_factory=list, init=False)
    translation_shift : int = field(init=False)
    area_correction : float = field(init=False)


@dataclass(eq=False)
class CompoundData(DadData, _CompoundDataDefaultsBase, _CompoundDataBase):
    """
    Parameter order: hplc_system_tag, experiment, gradient, wl_high_pass, wl_low_pass, detector_limit
    """

    def __post_init__(self, gradient, wl_high_pass, wl_low_pass):
        if gradient is None:
            raise TypeError("__init__ missing 1 required argument: 'gradient'")
        super().__post_init__(wl_high_pass, wl_low_pass)
        self._subtract_baseline(gradient)

    def _trim_data(self, length):
        """Trims the data in the time dimension to the length provided"""
        self.data, self.time = trim_data(data=self.data, time=self.time, length=length)
    
    def _subtract_baseline(self, gradient):
        """Subtracts the baseline of the gradient numpy array from self.data."""
        self._trim_data(gradient.data.shape[1])
        self.data = self.data - gradient.data[:, :self.data.shape[1]]

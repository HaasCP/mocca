# flake8: noqa
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 13:16:51 2021

@author: haascp
"""


from dataclasses import dataclass, field
from typing import Dict, List
from mocca.peak.models import BasePeak

import numpy as np
from mocca.dtypes.process_funcs import trim_data, get_compound_names, get_peaks, absorbance_to_array
from mocca.dtypes.process_gradientdata import bsl_als
from mocca.chemstation_api import read_csv_agilent, tidy_csv_agilent
from mocca.labsolutions_api import read_txt_shimadzu
from mocca.read_file_data import read_csv, apply_filter


@dataclass
class DADData():
    # actual class inputs
    hplc_system_tag : str
    path : str
    idx : int

    # set during initialization
    detector_limit : int = field(init=False)
    data : np.ndarray = field(init=False)
    time : np.ndarray = field(init=False)
    wavelength : np.ndarray = field(init=False)

    def __post_init__(self):
        self._set_detector_limit()
        self._read_data()

    def _set_detector_limit(self):
        if self.hplc_system_tag == 'chemstation':
            self.detector_limit = 2000
        elif self.hplc_system_tag == 'testing':
            self.detector_limit = np.inf
        elif self.hplc_system_tag == 'labsolutions':
            self.detector_limit = 2000
        else:
            raise AttributeError("HPLC System Tag {} not supported!".format(self.hplc_system_tag))

    def _read_data(self, wl_high_pass=None):
        if self.hplc_system_tag == 'chemstation':
            df = read_csv_agilent(self.path)
            df = apply_filter(tidy_csv_agilent(df, wl_high_pass))
        elif self.hplc_system_tag == 'testing':
            df = read_csv(self.path, encoding='utf-8')
            df = apply_filter(tidy_csv_agilent(df))
        elif self.hplc_system_tag == 'labsolutions':
            df = read_txt_shimadzu(self.path, wl_high_pass)
            df = apply_filter(df)

        self.data = absorbance_to_array(df)
        self.time = df.time.unique()
        self.wavelength = df.wavelength.unique()


@dataclass
class GradientData(DADData):
    def __post_init__(self):
        super().__post_init__(self)
        self.data = bsl_als(self.data)


@dataclass
class CompoundData(DADData):
    # input dictionary for compounds mapping compound_id to concentration
    compound_input : Dict

    # data properties to be set by class, do not initialize
    peaks : List[BasePeak] = field(default_factory=list)
    warnings : List[str] = field(default_factory=list)
    translation_shift : int = None
    area_correction : float = None

    def trim_data(self, length):
        """Trims the data in the time dimension to the length provided"""
        self.data, self.time = trim_data(data=self.data, time=self.time, length=length)

    def get_compound_names(self):
        # returns all compounds in input with nonzero concentration
        return get_compound_names(self.compound_input)

    def get_compound_concentration(self, key):
        return self.compound_input[key]

    def subtract_baseline(self, gradient_data):
        """Subtracts the baseline of the gradient numpy array from self.data."""
        self.trim_data(gradient_data.shape[1])
        self.data = self.data - gradient_data[:, :self.data.shape[1]]

    def get_peaks(self, absorbance_threshold, peaks_high_pass=None, peaks_low_pass=None, expand_peaks=True):
        """Generates a list of peaks in the spectra, in dictionary format"""
        return get_peaks(self, absorbance_threshold, peaks_high_pass, peaks_low_pass, expand_peaks)

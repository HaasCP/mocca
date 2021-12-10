#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 13:16:51 2021

@author: haascp
"""

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np
from mocca.dad_data.utils import absorbance_to_array, apply_filter, trim_data
from mocca.dad_data.process_gradientdata import bsl_als

from mocca.dad_data.apis.chemstation_api import read_csv_agilent, tidy_df_agilent
from mocca.dad_data.apis.labsolutions_api import read_txt_shimadzu


@dataclass()
class DadData():
    # actual class inputs
    hplc_system_tag : str
    path : str

    # set during initialization
    detector_limit : int = field(init=False)
    data : np.ndarray = field(init=False)
    time : np.ndarray = field(init=False)
    wavelength : np.ndarray = field(init=False)

    def __post_init__(self):
        self._set_detector_limit()
        self._read_data()

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            # don't attempt to compare against unrelated types
            raise ValueError("Both DAD datasets must be of the same type!")
        return self.path == other.path

    def _set_detector_limit(self):
        if self.hplc_system_tag == 'chemstation':
            self.detector_limit = 2000
        elif self.hplc_system_tag == 'labsolutions':
            self.detector_limit = 2000
        else:
            raise AttributeError("HPLC System Tag {} not supported!".format(self.hplc_system_tag))

    def _read_data(self, wl_high_pass=None):
        if self.hplc_system_tag == 'chemstation':
            df = read_csv_agilent(self.path)
            df = apply_filter(tidy_df_agilent(df, wl_high_pass))
        elif self.hplc_system_tag == 'labsolutions':
            df = read_txt_shimadzu(self.path, wl_high_pass)
            df = apply_filter(df)

        self.data = absorbance_to_array(df)
        self.time = df.time.unique()
        self.wavelength = df.wavelength.unique()


@dataclass(eq=False)
class GradientData(DadData):
    def __post_init__(self):
        super().__post_init__()
        self.data = bsl_als(self.data)


@dataclass(eq=False)
class CompoundData(DadData):
    # input dictionary for compounds mapping compound_id to concentration
    compound_input : Dict
    gradient : GradientData

    # data properties to be set by class, do not initialize
    warnings : List[str] = field(default_factory=list)
    translation_shift : int = None
    area_correction : float = None
    
    def __post_init__(self):
        super().__post_init__()
        self._subtract_baseline()

    def _trim_data(self, length):
        """Trims the data in the time dimension to the length provided"""
        self.data, self.time = trim_data(data=self.data, time=self.time, length=length)
    
    def _subtract_baseline(self):
        """Subtracts the baseline of the gradient numpy array from self.data."""
        self._trim_data(self.gradient.data.shape[1])
        self.data = self.data - self.gradient.data[:, :self.data.shape[1]]

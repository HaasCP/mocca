#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 14:57:25 2021

@author: haascp
"""
from dataclasses import dataclass
from typing import Optional


@dataclass()
class Settings():
    """
    Data container to store all user given data analysis settings.
    """
    hplc_system_tag : str
    detector_limit : Optional[float] = None
    absorbance_threshold : Optional[float] = 500
    wl_high_pass : Optional[float] = None
    wl_low_pass : Optional[float] = None
    peaks_high_pass : Optional[float] = None
    peaks_low_pass : Optional[float] = None
    spectrum_correl_thresh : Optional[float] = 0.95  # Value between 0 and 1
    relative_distance_thresh : Optional[float] = 0.01  # Value between 0 and 1

    def __post_init__(self):
        if self.detector_limit is None:
            if self.hplc_system_tag == 'chemstation':
                self.detector_limit = 2000
            elif self.hplc_system_tag == 'labsolutions':
                self.detector_limit = 2000
            elif self.hplc_system_tag == 'empower':
                self.detector_limit = 2000
            elif self.hplc_system_tag == 'custom':
                self.detector_limit = float("inf")
            else:
                raise AttributeError(f"HPLC System Tag {self.hplc_system_tag} n"
                                     "ot supported!")
        if self.hplc_system_tag == 'custom' and (self.wl_high_pass or
                                                 self.wl_low_pass):
            raise AttributeError("Wavelength high and low pass filters are not "
                                 "supported for custom data. Provide already "
                                 "trimmed data!")

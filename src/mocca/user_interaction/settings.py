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
            else:
                raise AttributeError("HPLC System Tag {} not supported!".format(self.hplc_system_tag))
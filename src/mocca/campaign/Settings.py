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
    detector_limit : Optional[float] = None
    absorbance_threshold : Optional[float] = None
    wl_high_pass : Optional[float] = None
    wl_low_pass : Optional[float] = None
    peaks_high_pass : Optional[float] = None
    peaks_low_pass : Optional[float] = None
    spectrum_correl_thresh : Optional[float] = None  # Value between 0 and 1
    relative_distance_thresh : Optional[float] = None  # Value between 0 and 1
    
    def update(self, detector_limit, absorbance_threshold, wl_high_pass,
               wl_low_pass, peaks_high_pass, peaks_low_pass,
               spectrum_correl_thresh, relative_distance_thresh):
        self.detector_limit = detector_limit
        self.absorbance_threshold = absorbance_threshold
        self.wl_high_pass = wl_high_pass
        self.wl_low_pass = wl_low_pass
        self.peaks_high_pass = peaks_high_pass
        self.peaks_low_pass = peaks_low_pass
        self.spectrum_correl_thresh = spectrum_correl_thresh
        self.relative_distance_thresh = relative_distance_thresh

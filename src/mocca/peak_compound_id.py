#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 12:06:57 2021

@author: haascp
"""

import mocca.peak
import numpy as np


class CompoundIdPredictor():
    def _average_peak_spectrum(self, peak):
        """
        Calculates mean spectrum over peak from left to right border.
        """
        return np.average(peak.dataset.data[:, peak.left:peak.right+1],
                          axis=0).tolist()
        
    def predict_compound_id(self, peak, component_db, len_time_vec):
        peak_spectrum = self._average_peak_spectrum(peak)
        for component in component_db:
            similarity = np.corrcoef(component.spectrum,
                                     peak_spectrum)[1, 0]
        
        
        
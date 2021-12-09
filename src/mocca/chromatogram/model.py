# flake8: noqa
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 13:22:39 2021

@author: haascp
"""
from typing import List

from mocca.peak.models import PreprocessedPeak
from mocca.chromatogram.utils import check_same_dataset

class Chromatogram():
    def __init__(self, peaks: List[PreprocessedPeak]):
        for peak in peaks:
            check_same_dataset(peak, peaks[0])
        self.peaks = peaks
    
class InitializationChromatogram(Chromatogram):
    def __init__(self, istd_key=None):
        self.istd_key = istd_key

    def process_peaks(self):
        compounds = get_compounds(dataset)
        if not any(peak.matches['compound_id'].startswith(compounds[0]) for peak in self.peaks):
            pass
        pass

class ReactionChromatogram(Chromatogram):
    def assign_peaks(self):
        max_similarity_peak = self.peaks[0]
        for peak in self.peaks:
            sorted(peak.matches, reverse=True,
                                   key=lambda dic: dic['spectrum_correl_coef'])
            if peak.matches > max_similarity_peak:
                pass
    
    def process_peaks():
        pass

# flake8: noqa
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 13:22:39 2021

@author: haascp
"""

from mocca.chromatogram.utils import check_same_dataset

class Chromatogram():
    def __init__(self):
        self.peaks = []
        self.dataset = None
    
    def insert_peak(self, peak):
        if self.dataset:
            check_same_dataset(peak, self.peaks[0])
        else:
            self.dataset = peak.dataset
        self.peaks.append(peak)

    def __iter__(self):
        """Yields all items inside the database.
           Allows statements such as "for component in component_database:" """
        for peak in self.peaks:
            yield peak
    
class InitializationChromatogram(Chromatogram):
    def __init__(self, istd_key=None):
        self.istd_key = istd_key

    def __contains__(self, compound_id):
        """Allows statements such as "'Component 1' in component_database"
        to see if that Component is inside the database"""
        return compound_id in [item.compound_id for item in self.items]

    def __getitem__(self, compound_id):
        """Allows statements such as "component_database['Component 1']" to access
           that Component"""
        if compound_id in self:
            return next(item for item in self.items if
                        item.compound_id == compound_id)
        else:
            raise AttributeError("{} not found in Database!"
                                 "".format(compound_id))

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

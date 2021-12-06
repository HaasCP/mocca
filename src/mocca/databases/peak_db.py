#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 08:28:12 2021

@author: haascp
"""

from mocca.peak.models import ProcessedPeak

import logging


class PeakDatabase():
    def __init__(self):
        self.peaks = []
    
    def __iter__(self):
        """Yields all items inside the database.
           Allows statements such as "for component in component_database:" """
        for peak in self.peaks:
            yield peak

    def __contains__(self, peak):
        """Allows statements such as "'Component 1' in component_database"
        to see if that Component is inside the database"""
        return True in [peak == db_peak for db_peak in self.peaks]

    def insert_peak(self, new_peak):
        if type(new_peak) != ProcessedPeak:
            raise TypeError("Warning: Data {} is not of the mocca ProcessedPeak type."
                            "Given data not inserted in the database".format(new_peak))
        elif new_peak in self:
            for i, peak in enumerate(self.peaks):
                if peak == new_peak:
                    self.peaks[i] = new_peak
                logging.warning("Warning: Peak \n {} \n already exists in database. "
                                "New peak replaces the "
                                "old one.".format(new_peak))
        else:
            self.peaks.append(new_peak)
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 08:28:12 2021

@author: haascp
"""

from mocca.peak.models import ProcessedPeak

import logging
from typing import Optional, List


class PeakDatabase():
    """
    Database class to store and organize peaks of HPLC-DAD data.
    """
    def __init__(self, peaks: Optional[List[ProcessedPeak]] = None):
        """
        Takes a list of peaks as an optional argument. If peaks are passed,
        they are checked for the processed peak datatype. If no list of peaks
        is passed, an empty list is instanstiated.
        """
        if peaks:
            if any(type(peak) != ProcessedPeak for peak in peaks):
                raise TypeError("Warning: List contains data which is not of the"
                                " mocca ProcessedPeak type."
                                "Given data not inserted in the database.")
            else:
                self.peaks = peaks
        else:
            self.peaks = []

        self.unknown_counter = 0

    def __iter__(self):
        """Yields all items inside the database.
           Allows statements such as "for component in component_database:" """
        for peak in self.peaks:
            yield peak

    def __contains__(self, peak):
        """Allows statements such as "peak in peak_db"
        to see if that Peak is inside the database"""
        return True in [peak == db_peak for db_peak in self.peaks]

    def update_unknown_counter(self):
        """
        Updates the unknown counter which gives the highest trailing number in
        the peak database of a peak assigned with a compound_id starting with
        'unknown_'.
        """
        cur_count = 0
        for peak in self:
            if peak.compound_id.startswith("unknown_"):
                num = int(peak.compound_id[8:])
                if num > cur_count:
                    cur_count = num
        self.unknown_counter = cur_count

    def increment_unknown_counter(self):
        """
        Increments the unknown counter by one.
        """
        self.unknown_counter += 1

    def insert_peak(self, new_peak):
        """
        Inserts a peak in the database given it is of the processed peak type.
        If a peak with the same borders in the same dataset already exists,
        it will be overwritten.
        """
        if type(new_peak) != ProcessedPeak:
            raise TypeError("Warning: Data {} is not of the mocca ProcessedPeak type."
                            "Given data not inserted in the database".format(new_peak))
        elif new_peak in self:
            for i, peak in enumerate(self.peaks):
                if peak == new_peak:
                    self.peaks[i] = new_peak
                logging.warning("Warning: Peak with maximum at {} already exists "
                                "in database. New peak replaces the "
                                "old one.".format(new_peak.maximum))
        else:
            self.peaks.append(new_peak)
        # update the unknown counter
        self.update_unknown_counter()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 08:28:12 2021

@author: haascp
"""

from mocca.peak import Peak
import logging


class PeakDatabase():
    def __init__(self):
        self.peaks = []

    def insert_peak(self, new_peak):
        if type(new_peak) != Peak:
            logging.warning("Warning: Data {} is not of the mocca Peak type."
                            "Given data not inserted in the database".format(new_peak))
            return
        else:
            for peak in self.peaks:
                if new_peak == peak:
                    logging.warning("Warning: Peak \n {} \n already exists in database "
                                    "as Peak \n {} \n New peak replaces the "
                                    "old one.".format(new_peak, peak))
                    self.peaks = [p for p in self.peaks if not p == peak]
            self.peaks.append(new_peak)
            logging.debug("Peak {} added to the database.".format(new_peak))

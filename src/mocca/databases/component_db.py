#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:32:36 2021

@author: haascp
"""
from mocca.databases.base import BaseDatabase
from mocca.databases.utils import get_filtered_peaks, get_filtered_peaks_by_compound
from mocca.component.funcs import create_component_from_peaks


class ComponentDatabase(BaseDatabase):
    def __init__(self):
        self.unknown_counter = 0

    def increment_unkown_counter(self):
        self.unknown_counter += 1

    def update_unknown_counter(self, peak_database):
        cur_count = 0
        for peak in peak_database:
            if peak.compound_id.startswith("unknown_"):
                num = int(peak.compound_id[8:])
                if num > cur_count:
                    cur_count = num
        return cur_count

    def update(self, peak_database, peak_filter_function=None):
        """
        Creates components from the given peak database. Optionally, a condition
        can be given to filter peaks.
        """
        # clear database to fill it with components
        self.delete_all_items()
        
        compound_dict = get_filtered_peaks_by_compound(peak_database.peaks, 
                                                       peak_filter_function)

        # create components out of compound dict
        for compound_id, peaks in compound_dict.items():
            component = create_component_from_peaks(peaks)
            self.insert_item(component)
        
        # update the unknown counter
        self.update_unknown_counter(peak_database)

    # for addition of unknown compounds in reaction runs if initialization runs
    # are not available anymore
    def insert_by_compound_id(self, peak_database, compound_id, 
                              peak_filter_function=None):
        """
        Inserts component in existing component list. If component with given
        compound_id already exists, it will be overwritten.
        """
        if compound_id in self:
            self.items = [item for item in self.items if
                          item.compound_id is not compound_id]
        
        filtered_peaks = get_filtered_peaks()
        
        compound_peaks = [peak for peak in filtered_peaks if
                          peak.compound_id == compound_id]

        component = create_component_from_peaks(compound_peaks)
        self.insert_item(component)

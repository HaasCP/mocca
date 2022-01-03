#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:32:36 2021

@author: haascp
"""
from mocca.components.utils import get_filtered_peaks, get_filtered_peaks_by_compound
from mocca.components.quali_funcs import create_quali_component


class BaseDatabase():
    """
    Base class for component databases with the unique constrained primary key
    compound_id.
    """
    def __init__(self):
        """
        Instantiates an empty list as the data container of the database
        """
        self.items = []

    def __iter__(self):
        """Yields all items inside the database.
           Allows statements such as "for component in component_database:" """
        for item in self.items:
            yield item

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

    def insert_item(self, item, unique=True):
        """
        Inserts a new item to the database. If the unique argument is True,
        it checks if the item already exists in the databse.
        """
        if unique:
            if item in self.items:
                raise NotImplementedError
            else:
                self.items.append(item)
        else:
            self.items.append(item)

    def delete_all_items(self):
        """
        Clears all components out of database.
        """
        self.items = []


class QualiComponentDatabase(BaseDatabase):

    def update(self, peak_database, peak_filter_function=None):
        """
        Creates components from the given peak database. Optionally, a condition
        can be given to filter peaks.
        """
        # clear database to fill it with components
        self.delete_all_items()
        compound_dict = get_filtered_peaks_by_compound(peak_database,
                                                       peak_filter_function)

        # create components out of compound dict
        for compound_id, peaks in compound_dict.items():
            component = create_quali_component(peaks)
            self.insert_item(component)

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
                          item.compound_id != compound_id]

        filtered_peaks = get_filtered_peaks(peak_database, peak_filter_function)

        compound_peaks = [peak for peak in filtered_peaks if
                          peak.compound_id == compound_id]
        if compound_peaks:
            component = create_quali_component(compound_peaks)
            self.insert_item(component)

class QuantComponentDatabase(BaseDatabase):
    def update(self, peak_database, peak_filter_function=None):
        pass

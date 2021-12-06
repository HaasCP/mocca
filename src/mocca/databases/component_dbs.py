#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:32:36 2021

@author: haascp
"""
import numpy as np
import logging

from mocca.databases.utils import get_filtered_peaks, get_filtered_peaks_by_compound
from mocca.databases.component_funcs import create_quali_component


class BaseDatabase():
    """
    Base class for component databases with the unique constrained primary key 
    compound_id.
    """
    def __init__(self):
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
        self.unknown_counter = cur_count

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
                          item.compound_id != compound_id]
        
        filtered_peaks = get_filtered_peaks(peak_database, peak_filter_function)
        
        compound_peaks = [peak for peak in filtered_peaks if
                          peak.compound_id == compound_id]
        if compound_peaks:
            component = create_quali_component(compound_peaks)
            self.insert_item(component)

class QuantiComponentDatabase():
    def __init__(self):
        self.database = {}

    def add_compound_concentration(self, compound_id, conc_info):
        """
        Adds tuple conc_info = (integral, concentration) to the component
        database under compound name name. Creates a new dictionary entry
        if that compound was originally not present.
        """
        if compound_id not in self:
            self.database[compound_id] = []
        self.database[compound_id].append((conc_info[0], conc_info[1]))

    def __contains__(self, compound_id):
        """Allows statements such as "'A' in quantification_database"
        to see if that Component is inside the quantification database"""
        return compound_id in self.database

    def quantify_peak(self, integral, compound_id):
        """
        Quantifies specified peak, given peak integral and compound id.
        """
        if compound_id not in self:
            raise Exception("Compound {} not found in Quantification Database!".format(compound_id))
        integrals, concentrations = zip(*self.database[compound_id])
        return np.linalg.lstsq(np.reshape(integrals, (-1, 1)), concentrations, rcond=0)[0][0] * integral
    
    def quantify_peak(self, quantification_database):
        """
        Computes the concentration of the compound in this peak.

        The attributes self.integral and self.compound_id must be set beforehand
        (through functions self.integrate_peak() and self.check_database())
        in order to quantify.

        Parameters
        ----------
        quantification_database : micdrop.QuantificationDatabase
            A quantification database containing all seen components.

        Raises Exception if the attributes self.integral or self.compound_id are
        not set. Prints a text warning if self.pure is not set.

        Modifies
        --------
        self.concentration : sets concentration to that predicted by integral
        """
        if self.integral is None or self.compound_id is None:
            raise Exception("Must set peak integral and compound_id before \
                            attempting to quantify!")

        if self.pure is None or not self.pure:
            logging.warning("Warning: Running quantify_peak() on impure peak \
                            {}.".format(self.idx))

        self.concentration = quantification_database.quantify_peak(self.integral,
                                                                   self.compound_id)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 09:05:19 2021

@author: haascp
"""
import logging
import dill

# init of campign instance
from mocca.peak.database import PeakDatabase
from mocca.components.databases import QualiComponentDatabase
from mocca.components.databases import QuantComponentDatabase

# process functions
from mocca.campaign.process_funcs import process_compound_experiments
from mocca.campaign.process_funcs import process_experiments
from mocca.campaign.process_funcs import process_gradients

# reporting
from mocca.report.main import report


class HplcDadCampaign():
    """
    Main parent class for HPLC-DAD experimental campaigns. We expect the
    gradient to stay constant over the campaign.
    """
    def __init__(self, autosave_path=None):
        self.autosave_path = autosave_path
        self.hplc_inputs = []
        self.settings = None
        self.peak_db = PeakDatabase()
        self.quali_comp_db = QualiComponentDatabase()
        self.quant_comp_db = QuantComponentDatabase()
        self.chroms = []
        self.warnings = []

    def _reset_campaign(self):
        """
        Resets the campaign to a state so that a full campaign data processing
        can take place.
        """
        self.peak_db = PeakDatabase()
        self.quali_comp_db = QualiComponentDatabase()
        self.quanti_comp_db = QuantComponentDatabase()
        self.compound_chroms = []
        self.chroms = []
        self.bad_chroms = []
        self.warnings = []
        for hplc_input in self.hplc_inputs:
            hplc_input.processed = False
            if hplc_input.gradient:
                hplc_input.gradient.dataset = None

    def save_campaign(self, path='hplc_dad_campaign.pkl', remove_raw_data=False):
        """
        Saves campaign object as pkl file. If remove_raw_data is True, all raw
        data are removed before saving to reduce file sizes.
        """
        if remove_raw_data:
            for chrom in self.chroms:
                chrom.dataset.time = []
                chrom.dataset.wavelength = []
                chrom.dataset.data = []
            for hplc_input in self.hplc_inputs:
                hplc_input.gradient.dataset.time = []
                hplc_input.gradient.dataset.wavelength = []
                hplc_input.gradient.dataset.data = []
        with open(path, 'wb') as file:
            dill.dump(self.__dict__, file)

    def load_campaign(self, path='hplc_dad_campaign.pkl'):
        """
        Loads campaign object which was saved as pkl file.
        """
        with open(path, 'rb') as file:
            self.__dict__.update(dill.load(file))

    def add_hplc_input(self, hplc_input):
        """
        All reaction concentration are metadata and should be trated like
        reaction temperature etc. Only give compound_conc if standard.
        If compound_id given: Add new peak to component db and update it
        If conc given: Add peak to quanti_component with same compound_id and update
        Store user input concs as negative conc in peak
        """
        for i, exp in enumerate(self.hplc_inputs):
            if exp.path == hplc_input.path:
                del self.hplc_inputs[i]
                logging.warning("CampaignWarning: Path {} was already used for "
                                "different HPLC input. New HPLC input replaces "
                                "the old one.".format(hplc_input.path))
        self.hplc_inputs.append(hplc_input)

        if self.autosave_path:
            self.save_campaign(path=self.autosave_path)

    def process_all_hplc_input(self, settings):
        """
        This function sets all expeirments of the cmapaign to the unprocessed
        state and processes all given hplc input. NOTE:
        This function has to be run if a new compound is added to the component
        database via compound experiment so that all peaks are assigned consistently
        """
        self.settings = settings
        self._reset_campaign()

        process_gradients(self.hplc_inputs, self.settings)

        chroms = process_compound_experiments(
            self.hplc_inputs,
            self.peak_db,
            self.quali_comp_db,
            self.quant_comp_db,
            self.settings
            )
        self.chroms.extend(chroms)

        chroms = process_experiments(
            self.hplc_inputs,
            self.peak_db,
            self.quali_comp_db,
            self.quant_comp_db,
            self.settings
            )
        self.chroms.extend(chroms)

    def process_new_hplc_input(self):
        """
        Only unprocessed runs are analyzed. No compound runs are allowed. Settings
        can only be changed via provess_all_hplc_input.
        """
        process_gradients(self.hplc_inputs, self.settings)

        chroms = process_experiments(
            self.hplc_inputs,
            self.peak_db,
            self.quali_comp_db,
            self.quant_comp_db,
            self.settings
            )
        self.chroms.extend(chroms)

    def generate_reports(self, path):
        """
        Consolidates all report functions in one function.
        """
        report(self, path)

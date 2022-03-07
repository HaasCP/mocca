#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 09:05:19 2021

@author: haascp
"""
import logging

# init of campign instance
from mocca.peak.database import PeakDatabase
from mocca.components.databases import QualiComponentDatabase
from mocca.components.databases import QuantComponentDatabase

# process functions
from mocca.campaign.process_funcs import process_compound_experiments
from mocca.campaign.process_funcs import process_experiments
from mocca.campaign.process_funcs import process_gradients
from mocca.campaign.utils import save_instance


class HplcDadCampaign():
    """
    Main parent class for HPLC-DAD experimental campaigns. We expect the
    gradient to stay constant over the campaign.
    """
    def __init__(self, autosave_path=None):
        self.autosave_path = autosave_path
        self.hplc_runs = []
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
        for hplc_input in self.hplc_runs:
            hplc_input.processed = False
            hplc_input.gradient.dataset = None

    def add_hplc_input(self, hplc_input):
        """
        All reaction concentration are metadata and should be trated like
        reaction temperature etc. Only give compound_conc if standard.
        If compound_id given: Add new peak to component db and update it
        If conc given: Add peak to quanti_component with same compound_id and update
        Store user input concs as negative conc in peak
        """
        for i, exp in enumerate(self.hplc_runs):
            if exp.path == hplc_input.path:
                del self.hplc_runs[i]
                logging.warning("CampaignWarning: Path {} was already used for "
                                "different HPLC input. New HPLC input replaces "
                                "the old one.".format(hplc_input.path))
        self.hplc_runs.append(hplc_input)

        if self.autosave_path:
            save_instance(self, self.autosave_path)

    def process_all_hplc_input(self, settings):
        """
        This function has to be run if a new compound is added to the component
        database via compound experiment so that all peaks are assigned consistently
        """
        self.settings = settings
        self._reset_campaign()

        process_gradients(self.hplc_runs, self.settings)

        chroms = process_compound_experiments(
            self.hplc_runs,
            self.peak_db,
            self.quali_comp_db,
            self.quant_comp_db,
            self.settings
            )
        self.chroms.extend(chroms)

        chroms = process_experiments(
            self.hplc_runs,
            self.peak_db,
            self.quali_comp_db,
            self.quant_comp_db,
            self.settings
            )
        self.chroms.extend(chroms)

# TODO : process_unprocessed_experiments

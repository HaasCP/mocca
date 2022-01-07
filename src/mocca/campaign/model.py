#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 09:05:19 2021

@author: haascp
"""
import logging

from mocca.peak.database import PeakDatabase
from mocca.components.databases import QualiComponentDatabase, QuantComponentDatabase

from mocca.campaign.settings import Settings
from mocca.campaign.process_funcs import (get_gradient,
                                          process_compound_experiments,
                                          process_experiments)
from mocca.campaign.utils import save_instance

# TODO: Go through files and clean up uncommented functions + documentation


class HplcDadCampaign():
    """
    Main parent class for HPLC-DAD experimental campaigns. We expect the gradient to stay
    constant over the campaign.
    """
    def __init__(self, hplc_system_tag, autosave_path=None):
        self.gradient = None
        self.autosave_path = autosave_path
        self.experiments = []
        self.settings = Settings(hplc_system_tag)
        self.peak_db = PeakDatabase()
        self.quali_comp_db = QualiComponentDatabase()
        self.quant_comp_db = QuantComponentDatabase()
        self.chroms = []
        self.warnings = []

    def _reset_campaign(self):
        self.peak_db = PeakDatabase()
        self.quali_comp_db = QualiComponentDatabase()
        self.quanti_comp_db = QuantComponentDatabase()
        self.compound_chroms = []
        self.chroms = []
        self.bad_chroms = []
        self.warnings = []
        for experiment in self.experiments:
            experiment.processed = False

        self.gradient = get_gradient(self.experiments, self.settings)

    def add_experiment(self, experiment):
        """
        All reaction concentration are metadata and should be trated like
        reaction temperature etc. Only give compound_conc if standard.
        If compound_id given: Add new peak to component db and update it
        If conc given: Add peak to quanti_component with same compound_id and update
        Store user input concs as negative conc in peak
        """
        if self.autosave_path:
            save_instance(self, self.autosave_path)

        for i, exp in enumerate(self.experiments):
            if exp.path == experiment.path:
                del self.experiments[i]
                logging.warning("CampaignWarning: Path {} was already used for "
                                "different experiment. New experiment replaces "
                                "the old one.".format(experiment.path))
        self.experiments.append(experiment)

    def process_all_experiments(self, detector_limit=None,
                                absorbance_threshold=500, wl_high_pass=None, 
                                wl_low_pass=None, peaks_high_pass=None, 
                                peaks_low_pass=None, spectrum_correl_thresh=0.9, 
                                relative_distance_thresh=0.01):
        """
        This function has to be run if a new compound is added to the component
        database via compound experiment so that all peaks are assigned consistently
        """
        self.settings.update(detector_limit, absorbance_threshold,
                             wl_high_pass, wl_low_pass, peaks_high_pass,
                             peaks_low_pass, spectrum_correl_thresh, 
                             relative_distance_thresh)
        self._reset_campaign()
        
        chroms = process_compound_experiments(
            self.experiments,
            self.gradient,
            self.peak_db,
            self.quali_comp_db,
            self.quant_comp_db,
            self.settings
            )
        self.chroms.extend(chroms)

        chroms = process_experiments(
            self.experiments,
            self.gradient,
            self.peak_db,
            self.quali_comp_db,
            self.quant_comp_db,
            self.settings
            )
        self.chroms.extend(chroms)
    
# TODO : process_unprocessed_experiments

     
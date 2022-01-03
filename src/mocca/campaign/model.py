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
                                          process_compound_experiments)
from mocca.campaign.utils import save_instance

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
        self.quanti_comp_db = QuantComponentDatabase()
        self.warnings = []

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

        for experiment in self.experiments:
            experiment.processed == False
        self.peak_db = PeakDatabase()
        self.quali_comp_db = QualiComponentDatabase()
        self.quanti_comp_db = QuantComponentDatabase()

        self.gradient = get_gradient(self.experiments, self.settings)
        
        compound_chroms = process_compound_experiments(self.experiments,
                                                       self.gradient,
                                                       self.peak_db,
                                                       self.quali_comp_db,
                                                       self.settings)
        return compound_chroms
"""   
        
        
        solvent_exps = [exp for exp in self.experiments if exp.compound.solvent]
        for exp in solvent_exps:
            chromatogram = process_solvent_exp(exp, self.hplc_system_tag, self.gradient,
                                               self.quali_comp_db, self.settings)
            
            
            
        # compound_exps =
        if new_exp.istd and new_exp.istd.key not in self.quali_comp_db:
            raise ValueError("Internal standard {} unknown in this campaign. "
                             "First add the internal standard as pure "
                             "compound in a separate run!".format(new_exp.istd.key))
        
        process_solvents()
        def order_experiments(self):
            solvent_exps = []
            # insert in peak_db, update quali_comp_db
            istd_exps = []
            # insert in peak_db, update quali_comp_db
            compound_exps = []
            # insert in peak_db, update quali_comp_db
            calibration_exps = []
            # assign impurities only now after all initializaiton is done
            # insert in peak_db, update quali_comp_db, update quanti_comp_db
            analysis_exps = []

        compound_data = CompoundData(self.hplc_system_tag, path, 
                                     self.gradient, {})
        chromatogram = pick_peaks
    
    def process_new_experiment(self, path, compound_id=None, solvent=False, istd=False, compound_conc=None,
                       istd_id=None, istd_conc=None):
        "
        Here, thresholds of the campaign are taken --> we don't have to process all
        "
        exps_to_process = [exp for exp in self.experiments if not exp.processed]
        if any (exp.compound for exp in exps_to_process):
            raise ValueError("If compound should be added to existing campaign, "
                             "the user has to run process_all_experiments function "
                             "to have consistent peak labelling over all experiments.")
        if not bad_data:
            update_all_databases(peak, component, compound)
        pass
    


    

    #TODO def load_peak_database()

    def preprocess_campaign(self, absorbance_threshold, wl_high_pass=None, 
                        wl_low_pass=None, peaks_high_pass=None, 
                        peaks_low_pass=None, spectrum_correl_thresh=0.9, 
                        relative_distance_thresh=0.01):
        for experiment in self.experiments:
            compound_data = CompoundData(self.hplc_system_tag, experiment['path'], 
                                         self.gradient, wl_high_pass, wl_low_pass, 
                                         experiment['compound_input'])
            chromatogram = pick_peaks(compound_data, absorbance_threshold, 
                                      peaks_high_pass, peaks_low_pass)
            
            1. expand
            2. check
            3. integrate
            4. correct
            5. resolve_impure (only impure peaks with component entry close by)
            6. match
            
            chromatogram = preprocess_chromatogram(chromatogram, self.istd_key,
                                                   self.quali_component_db, 
                                                   absorbance_threshold, 
                                                   self.detector_limit, 
                                                   spectrum_correl_thresh,
                                                   relative_distance_thresh)

            self.chromatograms.append(chromatogram)

# TODO develop compound_id assignment algorithms in each child class
"""

     
     
     
     
     
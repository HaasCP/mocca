#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 09:05:19 2021

@author: haascp
"""
import logging

from mocca.peak.database import PeakDatabase
from mocca.components.databases import QualiComponentDatabase, QuantComponentDatabase

from mocca.dad_data.models import GradientData, CompoundData
from mocca.dad_data.process_funcs import pick_peaks

from mocca.chromatogram.preprocessor import preprocess_chromatogram
from mocca.chromatogram.assign import add_quali_component

from mocca.campaign.experiment import Experiment
from mocca.campaign.settings import Settings
from mocca.campaign.utils import save_instance

class HplcDadCampaign():
    """
    Main parent class for HPLC-DAD experimental campaigns. We expect the gradient to stay
    constant over the campaign.
    """
    def __init__(self, hplc_system_tag, gradient_path, autosave_path=None):
        self.hplc_system_tag = hplc_system_tag
        self.gradient = self._create_gradient(gradient_path)
        self.autosave_path = autosave_path
        self.experiments = []
        self.settings = Settings()
        self.peak_db = PeakDatabase()
        self.quali_comp_db = QualiComponentDatabase()
        self.quanti_comp_db = QuantComponentDatabase()

    def _create_gradient(self, gradient_path):
        """Adds gradient data to the HPLC instance. Gradient data are required
        for all following data analysis procedures
        
        Parameters
        ----------
        path : str
            Path where the HPLC-DAD data are stored
        """
        #extra since unrelated to the others, no injection at all!
        return GradientData(self.hplc_system_tag, gradient_path)

    def add_experiment(self, path, istd, compound):
        """
        All reaction concentration are metadata and should be trated like
        reaction temperature etc. Only give compound_conc if standard.
        If compound_id given: Add new peak to component db and update it
        If conc given: Add peak to quanti_component with same compound_id and update
        Store user input concs as negative conc in peak
        """
        if self.autosave_path:
            save_instance(self, self.autosave_path)

        new_exp = Experiment(path, compound, istd)

        for i, exp in enumerate(self.experiments):
            if exp.path == new_exp.path:
                del self.experiments[i]
                logging.warning("CampaignWarning: Path {} was already used for "
                                "different experiment. New experiment replaces "
                                "the old one.".format(new_exp.path))

        self.experiments.append(new_exp)

    def process_all_experiments(self, absorbance_threshold=500, wl_high_pass=None, 
                                wl_low_pass=None, peaks_high_pass=None, 
                                peaks_low_pass=None, spectrum_correl_thresh=0.9, 
                                relative_distance_thresh=0.01):
        """
        Easiest case: only one peak
        More complicated cases: with istd and/or solvents
        Assignes highest unmatched peak with compound_id and puts it in component_db.
        All other given compound_ids must already be in component db and peaks must
        be found in chromatogram.
        --> function for one unknown peak!
        Order:
            1. solvent
            2. istd
            3. 
        """
        self.settings.update(self, absorbance_threshold, wl_high_pass, wl_low_pass,
                             peaks_high_pass, peaks_low_pass, spectrum_correl_thresh, 
                             relative_distance_thresh)
        
        solvent_exps = [exp for exp in self.experiments if exp.compound.solvent]
        for exp in solvent_exps:
            compound_data = CompoundData(self.hplc_system_tag, exp, self.gradient,
                                         self.settings.wl_high_pass, self.settings.wl_low_pass)
            chromatogram = pick_peaks(compound_data, absorbance_threshold, 
                                      peaks_high_pass, peaks_low_pass)
            chromatogram = preprocess_chromatogram(chromatogram, self.istd_key,
                                                   self.quali_component_db, 
                                                   absorbance_threshold, 
                                                   self.detector_limit, 
                                                   spectrum_correl_thresh,
                                                   relative_distance_thresh)
            
            
        
        
        process_solvents()
        def order_experiments(self):
            solvent_exps = []
            # insert in peak_db, update quali_comp_db
            istd_exps = []
            # insert in peak_db, update quali_comp_db
            compound_exps = []
            # insert in peak_db, update quali_comp_db
            calibration_exps = []
            # insert in peak_db, update quali_comp_db, update quanti_comp_db
            analysis_exps = []

        compound_data = CompoundData(self.hplc_system_tag, path, 
                                     self.gradient, {})
        chromatogram = pick_peaks
    
    def process_new_experiment(self, path, compound_id=None, solvent=False, istd=False, compound_conc=None,
                       istd_id=None, istd_conc=None):
        """
        Here, thresholds of the campaign are taken --> we don't have to process all
        """
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
            """
            1. expand
            2. check
            3. integrate
            4. correct
            5. resolve_impure (only impure peaks with component entry close by)
            6. match
            """
            chromatogram = preprocess_chromatogram(chromatogram, self.istd_key,
                                                   self.quali_component_db, 
                                                   absorbance_threshold, 
                                                   self.detector_limit, 
                                                   spectrum_correl_thresh,
                                                   relative_distance_thresh)

            self.chromatograms.append(chromatogram)

# TODO develop compound_id assignment algorithms in each child class

class CalibrationCampaign(HplcDadCampaign):
    
    # We can even get rid of analyte_key since they are given via experiment input
    
    def process_initialization(self, absorbance_threshold, wl_high_pass=None, 
                               wl_low_pass=None, peaks_high_pass=None, 
                               peaks_low_pass=None, spectrum_correl_thresh=0.9, 
                               relative_distance_thresh=0.01):
        """
        Easiest case: only one peak
        More complicated cases: with istd and/or solvents
        Assignes highest unmatched peak with compound_id and puts it in component_db.
        All other given compound_ids must already be in component db and peaks must
        be found in chromatogram.
        --> function for one unknown peak!
        Order:
            1. solvent
            2. istd
            3. 
        """
        self.thresholds = {
            'absorbance_threshold': absorbance_threshold
            }
        compound_data = CompoundData(self.hplc_system_tag, path, 
                                     self.gradient, {})
        chromatogram = pick_peaks

    def add_calibration_run():
        pass

    def process_peaks(self):
        #compounds = get_compounds(dataset)
        #if not any(peak.matches['compound_id'].startswith(compounds[0]) for peak in self.peaks):
        #    pass
        pass
    
    def _create_quantification_db(self, plots=True, title=True, warning_threshold=0.98):
        """
        Creates a list of components from AnalyteData, characterized by pairs 
        of (istd-corrected integral, concentration), for use in quantification.

        If a pure component has more than one peak, the one with the highest
        non-istd integral is used.

        After this pairing is made, a linear regression is performed and the 
        relevant coefficient is stored in the quantification_database dict.
        """
        
        db = {}
        for data in self.data_container:
            if self.istd_key:    
                istd_peaks = [peak for peak in data.peaks if 'compound_id' in peak and peak['compound_id'] == self.istd_key]
                if len(istd_peaks) != 1:
                    raise AttributeError("More or less than one ISTD peak found!")
                istd_peak = istd_peaks[0]
                if 'integral' in istd_peak:
                    data.area_correction = data.input[self.istd_key][0] / istd_peak['integral']
                compound_name = [compound_name for compound_name in data.get_compound_names() if compound_name != self.istd_key][0]
            else:
                compound_name = data.get_compound_names()[0]
            
            if compound_name not in db:
                db[compound_name] = []
            peak_integrals = [peak['integral'] for peak in data.peaks if peak.get('compound_id') == compound_name]
            if len(peak_integrals) == 0:
                print(f'Warning: no {compound_name} peak found in data; not adding to quantification database')
            else:
                db[compound_name].append((peak_integrals[0] * data.area_correction if hasattr(data, 'area_correction') else peak_integrals[0], 
                                      data.get_compound_concentration(compound_name)))
    
# class ReactionCampaign(HplcDadCampaign):
     # can overtake component dbs from calib campaign but starts own peak_db
     # component db is not updated but rather inserted by compound
     # Another possibility is to create component dbs from existing peak db
     
     
     


     
     
     
     
     
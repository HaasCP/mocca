#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 09:05:19 2021

@author: haascp
"""
import dill

from mocca.data_data.models import GradientData, CompoundData
from mocca.dad_data.process_funcs import pick_peaks
from mocca.chromatogram.funcs import preprocess_chromatogram

class HplcDadCampaign():
    """
    Main parent class for HPLC-DAD experimental campaigns. We expect the gradient to stay
    constant over the campaign.
    """
    def __init__(self, hplc_system_tag, gradient_path, istd_id=None):
        self.hplc_system_tag = hplc_system_tag
        self.istd_key = istd_id
        self.gradient = self._create_gradient(gradient_path)
        self.experiments = []
        self.chromatograms = []
        self.peak_database = None
        self.quali_component_db = None
        self.quanti_component_db = None

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

    def save_instance(self, path='data.pkl'):
        file = open(path, 'wb')
        for chromatogram in self.chromatograms:
            for peak in chromatogram.peaks:
                peak.dataset.data = []
        dill.dump(self.__dict__, file)
        file.close()
    
    def load_instance(self, path='data.pkl'):
        file = open(path, 'rb')
        self.__dict__.update(dill.load(file))
        file.close()

    def add_experiments(self, experiments):
        """
        compound_inputs: list of dicts with scalar values (no lists)
        """
        if not isinstance(experiments, list):
            raise TypeError("Compound inputs have to be given as a list of dictionaries.")
        if not all("path" in compound_input
                   for compound_input in experiments):
            raise AttributeError("For all added experiments, a path to the data "
                                 "has to be given under the dictionary key 'path'.")
        if self.istd_id and not all(self.istd_id in compound_input
                                    for compound_input in experiments):
            raise AttributeError("In campaigns with an internal standard, an internal "
                                 "standard input has to be given via the compound input dictionary.")
        self.experiments.extend(experiments)

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
            chromatogram = preprocess_chromatogram(chromatogram, 
                                                   self.quali_component_db, 
                                                   absorbance_threshold, 
                                                   self.detector_limit, 
                                                   spectrum_correl_thresh,
                                                   relative_distance_thresh)
            # TODO add PARAFAC routine
            self.chromatograms.append(chromatogram)

# TODO develop compound_id assignment algorithms in each child class

class InitializationCampaign(HplcDadCampaign):
    
    def __init__(self, analyte_ids=[]):
        #super_init
        self.analyte_keys = analyte_ids
    
    
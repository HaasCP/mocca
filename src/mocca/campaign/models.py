#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 09:05:19 2021

@author: haascp
"""
import dill

from mocca.dad_data.models import GradientData, CompoundData
from mocca.dad_data.process_funcs import pick_peaks
from mocca.chromatogram.preprocessor import preprocess_chromatogram

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

class InitializationCampaign(HplcDadCampaign):
    
    def __init__(self, analyte_ids=[]):
        #super_init
        self.analyte_keys = analyte_ids

    def process_peaks(self):
        compounds = get_compounds(dataset)
        if not any(peak.matches['compound_id'].startswith(compounds[0]) for peak in self.peaks):
            pass
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
    
    
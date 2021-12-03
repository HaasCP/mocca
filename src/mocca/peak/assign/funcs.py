#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 12:06:57 2021

@author: haascp
"""
from mocca.peak.models import AssignedPeak
from mocca.peak.utils import average_peak_spectrum
import numpy as np

def get_spectrum_correl_coef(peak, component):
    """
    """
    peak_spectrum = average_peak_spectrum(peak)
    return np.corrcoef(peak_spectrum, component.spectrum)[1, 0]

def get_relative_distance(peak, component):
    """
    """
    distance = abs(peak.maximum - component.maximum)
    return distance / len(peak.dataset.time)

def get_similarity_dict(peak, component_db):
    """
    """
    simil_by_comp = []
    for component in component_db:
        dic = {}
        dic['compound_id'] = component.compound_id
        dic['spectrum_correl_coef'] = get_spectrum_correl_coef(peak, component)
        dic['distance'] = abs(peak.maximum - component.maximum)
        dic['relative_distance'] = get_relative_distance(peak, component)
        simil_by_comp.append(dic)
    simil_by_comp = sorted(simil_by_comp, 
                           key=lambda dic: dic['spectrum_correl_coef'])

def predict_compound_id(peak, component_db, spectrum_correl_coef_thresh,
                        relative_distance_thresh, print_out):
    """
    """
    similarity_dict = get_similarity_dict(peak, component_db)
    if print_out:
        print(similarity_dict)
    matches = [d for d in similarity_dict if (d['spectrum_correl_coef'] >
                                              spectrum_correl_coef_thresh and
                                              d['relative_distance'] <
                                              spectrum_correl_coef_thresh)]
    if matches:
        matches = sorted(matches, key=lambda d: d['spectrum_correl_coef'])    
        return matches[0]['compound_id']
    else:
        return None

def get_compound_id(peak, component_db, spectrum_correl_coef_thresh,
                    relative_distance_thresh, print_out):
    """
    Checks if the peak is in the database, and if so, sets self.compound_id.
    See check_database().

    Modifies
    --------
    self.compound_id
        If a match is found, then the peak's compound_id attribute is set.
    """
    if peak.pure is True:
        predicted_compound_id = predict_compound_id(peak, component_db, 
                                                    spectrum_correl_coef_thresh,
                                                    relative_distance_thresh, 
                                                    print_out)
        if predicted_compound_id is None:
            component_db.increment_unkown_counter()
            return "unknown_" + str(component_db.unknown_counter)
        else:
            return predicted_compound_id
    else:
        return None

def assign_peak(checked_peak, component_db, spectrum_correl_coef_thresh,
                    relative_distance_thresh, print_compound_prediction):
    compound_id = get_compound_id(checked_peak, component_db,
                                  spectrum_correl_coef_thresh,
                                  relative_distance_thresh, 
                                  print_compound_prediction)

    return AssignedPeak(left = checked_peak.left,
                        right = checked_peak.right,
                        maximum = checked_peak.maximum,
                        dataset = checked_peak.dataset,
                        idx = checked_peak.idx,
                        saturation = checked_peak.saturation,
                        pure = checked_peak.saturation,
                        compound_id = compound_id)
    
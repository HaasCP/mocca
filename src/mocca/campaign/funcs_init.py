#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 09:49:07 2021

@author: haascp
"""
from operator import attrgetter

from mocca.dad_data.utils import get_compound_ids, get_compound_concentration
from mocca.peak.process import process_peak
    

def get_istd_chrom(chromatograms, istd_key):
    if istd_key:
        len_one_counter = 0
        for chrom in chromatograms:
            compound_ids = get_compound_ids(chrom.dataset)
            if len(compound_ids) == 1:
                len_one_counter += 1
                istd_chrom = chrom
        if len_one_counter == 0:
            raise Exception("The user has to provide a pure internal standard "
                            "experiment in order to use data analysis with "
                            "internal standard.")
        elif len_one_counter > 1:
            raise Exception("Only one experiment with one compound was expected. "
                            "{} experiments were found. For all "
                            "experiments containing other compounds than internal "
                            "standard, internal standard still has to be added in "
                            "order to allow for data analysis with internal "
                            "standard.".format(len_one_counter))
        else:
            return istd_chrom


def get_max_integral_peak(peaks):
    """
    Returns the peak with the maximum integral value in the given list of peaks.
    """
    if not all(hasattr(peak, 'integral') for peak in peaks):
        raise AttributeError("All given peaks must have integral attribute.")
    return max(peaks, key=attrgetter('integral'))


def process_residual_peaks(peaks, compound_id):
    """
    Assigns impurity compound ids to peaks.
    """
    peaks = sorted(peaks, key=lambda peak: peak.maximum)
    new_peaks = []
    impurity_counter = 0
    for peak in peaks:
        impurity_counter += 1
        new_peak = process_peak(peak, 
                                compound_id + "_impurity_" + str(impurity_counter))
        new_peaks.append(new_peak)
    return new_peaks


def get_processed_istd_chrom(chromatograms, istd_key):
    """
    Assigns istd_key to maximum integral peak and labels residual peaks as 
    impurities.
    """
    istd_chrom = get_istd_chrom(chromatograms, istd_key)
    istd_conc = get_compound_concentration(istd_chrom.dataset, istd_key)
    
    istd_peak = get_max_integral_peak(istd_chrom)
    residual_peaks = [peak for peak in istd_chrom if peak != istd_peak]

    istd_peak = process_peak(istd_peak, istd_key, istd_conc)
    residual_peaks = process_residual_peaks(residual_peaks, istd_key)
    
    istd_chrom.peaks = sorted([istd_peak] + residual_peaks,
                              key=lambda peak: peak.maximum)

    return istd_chrom




def get_max_conc_chroms_by_compound(chromatograms, istd_key):
    max_conc_runs = {}
    for chrom in chromatograms:
        compound_ids = get_compound_ids(chrom.dataset)
        for compound_id in compound_ids:
            if compound_id not in max_conc_runs:
                max_conc_runs[compound_id] = chrom
            elif (get_compound_concentration(chrom, compound_id) >
                  get_compound_concentration(max_conc_runs[compound_id].dataset, 
                                             compound_id)):
                max_conc_runs[compound_id] = chrom
    if istd_key:
        del max_conc_runs[istd_key]
    return max_conc_runs


def get_quali_db_init_chroms(chromatograms, istd_key):
    chroms = get_max_conc_chroms_by_compound(chromatograms, istd_key)
    if istd_key:
        istd_chrom = get_istd_chrom(chromatograms, istd_key)
        chroms[istd_key] = istd_chrom
    return chroms

        
                
        
        
        
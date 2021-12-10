# flake8: noqa
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 11:38:45 2021

@author: haascp
"""

"""
Init_campaign: 
    if input_compound_id not in component_db:
        peak_max_absorbance.compound_id = input_compound_id
    else:
        peak.get_compound_id
    component_db.update()

# react
processed_peaks = []
for peak in picked_peaks:
    peak = process_peak(PickedPeak)
    processed_peaks.append(peak)
remove_duplicate_compound_ids
for peak in processed_peaks:
    peak_db.insert_item(peak)
    if not update_component_db:
        if peak.compound_id not in component_db:
            component_db.insert_by_compound_id(peak_database, compound_id, 
                                      peak_filter_function=None)
if update_component_db:
    component_db.update(peak_database, compound_id, 
                              peak_filter_function=None)

# init
processed_peaks = []
for peak in picked_peaks:
    peak = process_peak(PickedPeak)
    processed_peaks.append(peak)
remove_duplicate_compound_ids
remove_unknowns(processed_peaks)
for peak in processed_peaks:
    peak_db.insert_item(peak)
component_db.update()

Impure peak always gets None in compound_id and concentration



Get list of chromatograms. Input can be accessed via peak.dataset

def get_preprocessed_chromatogram(self, component_db, absorbance_threshold=500,
                                  peaks_high_pass=None, peaks_low_pass=None,
                                  spectrum_correl_thresh = 0.9,
                                  relative_distance_thresh = 0.1):
    chromatogram = pick_peaks(self, absorbance_threshold, peaks_high_pass, 
                              peaks_low_pass)
    return preprocess_chromatogram(chromatogram, component_db,
                                   absorbance_threshold, 
                                   self.detector_limit, 
                                   spectrum_correl_thresh,
                                   relative_distance_thresh)



"""
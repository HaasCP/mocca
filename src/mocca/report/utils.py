#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 16:28:04 2022

@author: haascp
"""
import pandas as pd

def settings_to_df(settings):
    settings_dict = {
        'hplc_system_tag': settings.hplc_system_tag,
        'detector_limit': settings.detector_limit,
        'absorbance_threshold': settings.absorbance_threshold,
        'wl_high_pass': settings.wl_high_pass,
        'wl_low_pass': settings.wl_low_pass,
        'peaks_high_pass': settings.peaks_high_pass,
        'peaks_low_pass': settings.peaks_low_pass,
        'spectrum_correl_thresh': settings.spectrum_correl_thresh,
        'relative_distance_thresh': settings.relative_distance_thresh
        }
    settings_df = pd.DataFrame(settings_dict, index=[0])
    return settings_df
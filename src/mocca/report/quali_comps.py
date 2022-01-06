#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 17:29:16 2022

@author: haascp
"""
import os
import pandas as pd
import datapane as dp
from scipy.signal import find_peaks

import matplotlib.pyplot as plt 


def quali_comps_to_dict(comps):
    quali_comp_dict = {'compound_id': [],
                       'left': [],
                       'right': [],
                       'maximum': [],
                       'lambda_max': [],
                       'num_peaks': []}
                       
    for comp in comps:
        quali_comp_dict['compound_id'].append(comp.compound_id)
        quali_comp_dict['left'].append(comp.left)
        quali_comp_dict['right'].append(comp.right)
        quali_comp_dict['maximum'].append(comp.maximum)
        wls = comp.created_from[0].dataset.wavelength
        spectrum_maxima, _ = find_peaks(comp.spectrum)
        spectrum_maxima = [m for m in spectrum_maxima if comp.spectrum[m] > 1]
        lambda_max = [wls[i] for i in spectrum_maxima]
        quali_comp_dict['lambda_max'].append(lambda_max)
        quali_comp_dict['num_peaks'].append(len(comp.created_from))
    return quali_comp_dict

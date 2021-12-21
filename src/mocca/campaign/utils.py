#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:47:14 2021

@author: haascp
"""

import logging
import numpy as np


def suggest_initialization_runs(n_calib_dict, max_conc_dict=None, istd_key=None,
                                n_solvents=0):
    """Returns data frame of suggested runs for a standardized HPLC 
    initialization (initialization = calibration + database of analytes)
    
    Parameters
    ----------
    n_calib_df : pd.DataFrame, optional
        Number of calibration points for each analyte. Form:
        pd.DataFrame([[n(A1), n(A2), n(A3)]], 
                         columns=self.analyte_keys)

    max_conc_df : pd.DataFrame, optional
        Dataframe of the maximum absolute concentrations to be used in the
        calibration runs. Form:
            pd.DataFrame([[c(A1), c(A2), c(ISTD)]], 
                         columns=self.analyte_keys + [self.istd_key])
    
    n_solvents : int, optional
        Number of UV-Vis active solvents used in the reaction campaign. 
        For each solvent, a blank solvent run is added.
        
    export_path : str, optional
        Full path of the csv file where suggested runs are exported
    """

    if any(n_calib < 3 for n_calib in n_calib_dict.values()):
        logging.warning("Warning: Linear regression with less than three "
                        "calibration points not reliable. Increase n_calib_df "
                        "to >= 3 for quantification")

    calib_dicts = []
    calib_runs = []
    # add calibration runs for every analyte
    for analyte_id, n_runs in n_calib_dict.items():
         # create relative conentrations between 0 and 1 excluding 0
        conc_vec = np.linspace(1, 0, n_runs, endpoint = False).tolist()
        conc_vec = [round(num, 4) for num in conc_vec]
        conc_dict = {}
        conc_dict[analyte_id] = conc_vec
        calib_dicts.append(conc_dict)   
        calib_runs.append(n_runs)
    
    calib_dict = {}
    for k in set().union(*calib_dicts):
        calib_dict[k] = []
        for i, d in enumerate(calib_dicts):
            calib_dict[k].append(d.get(k, [0] * calib_runs[i]))
    
    for key in calib_dict.keys():
        calib_dict[key] = [item for sublist in calib_dict[key] for item in sublist]
    
    # add run with internal standard only
    if istd_key:
        for val in calib_dict.values():
            val.insert(0, 0)
        calib_dict[istd_key] = [1] * len(val)
    
    # add run with solvent only for UV-Vis active solvents
    if n_solvents > 0:
        for i in range(n_solvents):
            for val in calib_dict.values():
                val.insert(0, 0)
                
    # multiply relative concentrations with the given maximum concentrations
    if max_conc_dict:
        for analyte_id, max_conc in max_conc_dict.items():
            calib_dict[analyte_id] = [val * max_conc for val in calib_dict[analyte_id]]
    
    initialization_runs = calib_dict
    return initialization_runs
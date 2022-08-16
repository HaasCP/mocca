#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 09:10:24 2022

@author: haascp
"""

import logging
import numpy as np


def suggest_initialization_runs(n_calib_dict, max_conc_dict=None, istd_key=None,
                                n_solvents=0):
    """Returns data frame of suggested runs for a standardized HPLC
    initialization (initialization = calibration + database of analytes).
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
        conc_vec = np.linspace(1, 0, n_runs, endpoint=False).tolist()
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

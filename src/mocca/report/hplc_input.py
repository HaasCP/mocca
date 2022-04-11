#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 14:12:58 2022

@author: haascp
"""

import os
import pandas as pd
import datapane as dp


def exps_to_df(exps):
    """
    Transfers relevant information of HplcInput objects in a pandas df.
    """
    exp_dict = {'index': [],
                'file': [],
                'compound_key': [],
                'compound_conc': [],
                'compound_is_solvent': [],
                'compound_is_istd': [],
                'istd_keys': [],
                'istd_concs': [],
                'gradient_file': [],
                'processed': []}
    for i, exp in enumerate(exps):
        exp_dict['index'].append(i + 1)
        exp_dict['file'].append(os.path.basename(exp.path))
        if exp.compound:
            exp_dict['compound_key'].append(exp.compound.key)
            exp_dict['compound_conc'].append(exp.compound.conc)
            exp_dict['compound_is_solvent'].append(exp.compound.is_solvent)
            exp_dict['compound_is_istd'].append(exp.compound.is_istd)
        else:
            exp_dict['compound_key'].append(None)
            exp_dict['compound_conc'].append(None)
            exp_dict['compound_is_solvent'].append(None)
            exp_dict['compound_is_istd'].append(None)
        if exp.istd:
            istd_keys = [istd.key for istd in exp.istd]
            istd_concs = [istd.conc for istd in exp.istd]
            exp_dict['istd_keys'].append(istd_keys)
            exp_dict['istd_concs'].append(istd_concs)
        else:
            exp_dict['istd_keys'].append(None)
            exp_dict['istd_concs'].append(None)
        if exp.gradient:
            exp_dict['gradient_file'].append(os.path.basename(exp.gradient.path))
        else:
            exp_dict['gradient_file'].append(None)
        exp_dict['processed'].append(exp.processed)
    return pd.DataFrame(exp_dict)


def report_hplc_input(exps, report_path):
    """
    Main HPLC input report function.
    """
    if not exps:
        print("No HPLC input was given!")
        return
    this_dir, _ = os.path.split(__file__)
    mocca_icon_path = os.path.join(this_dir, "mocca_icon.png")
    exp_df = exps_to_df(exps)
    exp_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# HPLC input report"),
                dp.Media(file=mocca_icon_path),
                columns=2
            ),
            dp.Text("### Table: HPLC input as given by the user."),
            dp.DataTable(exp_df, label="experiment_table")
        ],
    )
    r = dp.Report(
        exp_page
    )
    r.save(path=os.path.join(report_path, "report_hplc_input.html"), open=True)

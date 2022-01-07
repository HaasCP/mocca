#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 14:12:58 2022

@author: haascp
"""

import os
import pandas as pd
import datapane as dp


def exps_to_dict(exps):
    exp_dict = {'file': [],
                'compound_key': [],
                'compound_conc': [],
                'compound_is_solvent': [],
                'compound_is_istd': [],
                'istd_keys': [],
                'istd_concs': [],
                'is_gradient': [],
                'processed': []}
    for exp in exps:
        exp_dict['file'].append(os.path.basename(exp.path))
        if exp.compound:
            exp_dict['compound_key'].append(exp.compound.key)
            exp_dict['compound_conc'].append(exp.compound.conc)
            exp_dict['compound_is_solvent'].append(exp.compound.solvent)
            exp_dict['compound_is_istd'].append(exp.compound.istd)
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
        exp_dict['is_gradient'].append(exp.gradient)
        exp_dict['processed'].append(exp.processed)
    return exp_dict


def exps_to_df(exps):
    exp_dict = exps_to_dict(exps)
    exp_df = pd.DataFrame(exp_dict)
    return exp_df


def report_experiments(exps, report_path):
    exp_df = exps_to_df(exps)
    exp_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# 1 Experiment report"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text("### Table: Experiments as given by the user."),
            dp.DataTable(exp_df, label="experiment_table")
        ],
    )
    r = dp.Report(
        exp_page
    )
    r.save(path=os.path.join(report_path, "report_experiments.html"), open=True)


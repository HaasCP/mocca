#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 09:36:47 2022

@author: haascp
"""
import os
import pandas as pd
import datapane as dp

from mocca.utils.visualization import plot_calibration_curves
from mocca.report.peaks import peaks_to_df


def get_max_score_version(comp):
    max_score = 0
    max_version = 'absolute'
    for version, score in comp.calib_scores.items():
        if score > max_score:
            max_score = score
            max_version = version
    return max_version

def quant_comps_to_dict(comps):
    quant_comp_dict = {'compound_id': [],
                       'version': [],
                       'calib_factor': [],
                       'R-squared': [],
                       'num_peaks': []}
                       
    for comp in comps:
        max_version = get_max_score_version(comp)
        quant_comp_dict['compound_id'].append(comp.compound_id)
        quant_comp_dict['version'].append(max_version)
        quant_comp_dict['calib_factor'].append(comp.calib_factors[max_version])
        quant_comp_dict['R-squared'].append(comp.calib_scores[max_version])
        quant_comp_dict['num_peaks'].append(len(comp.created_from))
    return quant_comp_dict


def comps_to_df(comps):
    comp_dict = quant_comps_to_dict(comps)
    comp_df = pd.DataFrame(comp_dict)
    return comp_df


def create_quant_comp_page(comp):
    calibration_curves = plot_calibration_curves(comp)
    plots = [dp.Plot(calib_curve) for calib_curve in calibration_curves]

    peaks_df = peaks_to_df(comp.created_from)
    return dp.Page(
        title=comp.compound_id,
        blocks=[
            dp.Group(
                dp.Text(f"## Details to component {comp.compound_id}"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text(f"### Figures: Calibration curves of component {comp.compound_id}.")] +
            plots +
           [ dp.Text("### Table: Peaks from which the component was created."),
            dp.DataTable(peaks_df, label=f"peak_table_{comp.compound_id}")
        ],        
    )


def report_quant_comps(quant_comp_db, report_path):
    comps = quant_comp_db.items
    comp_df = comps_to_df(comps)
    table_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# 4 Quantitative component database report"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text("### Table: Components in the quantitative component database."
                    "of the campaign."),
            dp.DataTable(comp_df, label="comp_table")
        ],        
    )
    component_pages = []
    for comp in quant_comp_db:
        page = create_quant_comp_page(comp)
        component_pages.append(page)
    r = dp.Report(
        table_page,
        *component_pages
    )
    r.save(path=os.path.join(report_path, "report_quant_comp_db.html"), open=True)


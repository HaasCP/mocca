#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 09:36:47 2022

@author: haascp
"""
import os
import pandas as pd
import datapane as dp

from mocca.visualization.calibration_plots import plot_calibration_curves
from mocca.report.peaks import peaks_to_df


def get_max_score_version(comp):
    """
    Returns the version of the calibration which has maximum R-squared.
    """
    max_score = 0
    max_version = 'absolute'
    for version, score in comp.calib_scores.items():
        if score > max_score:
            max_score = score
            max_version = version
    return max_version


def quant_comps_to_df(comps):
    """
    Transfers relevant information from quantitative components in a pandas df.
    """
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
    return pd.DataFrame(quant_comp_dict)


def create_quant_comp_page(comp):
    """
    Creates a report page for the given quantitative component.
    """
    calibration_curves = plot_calibration_curves(comp)

    peaks_df = peaks_to_df(comp.created_from)
    return dp.Page(
        title=comp.compound_id,
        blocks=[
            dp.Group(
                dp.Text(f"## Details to component {comp.compound_id}"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
                ),
            dp.Text("### Figures: Calibration curves of component "
                    f"{comp.compound_id}."),
            dp.Group(
                *calibration_curves,
                columns=2
                ),
            dp.Text("### Table: Peaks from which the component was created."),
            dp.DataTable(peaks_df, label=f"peak_table_{comp.compound_id}")
        ],
    )


def report_quant_comps(quant_comp_db, report_path):
    """
    Creates html report for the quantitative component database.
    """
    this_dir, _ = os.path.split(__file__)
    mocca_icon_path = os.path.join(this_dir, "mocca_icon.png")
    comps = quant_comp_db.items
    comp_df = quant_comps_to_df(comps)
    table_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# Quantitative component database report"),
                dp.Media(file=mocca_icon_path),
                columns=2
            ),
            dp.Text("### Table: Components in the quantitative component database "
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

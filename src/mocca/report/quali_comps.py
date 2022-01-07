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

from mocca.utils.visualization import plot_1D_data
from mocca.report.peaks import peaks_to_df


def quali_comps_to_dict(comps):
    quali_comp_dict = {'compound_id': [],
                       'left': [],
                       'right': [],
                       'maximum': [],
                       'lambda_max': [],
                       'num_peaks': []}
                       
    for comp in comps:
        times = comp.created_from[0].dataset.time
        quali_comp_dict['compound_id'].append(comp.compound_id)
        quali_comp_dict['left'].append(times[comp.left])
        quali_comp_dict['right'].append(times[comp.right])
        quali_comp_dict['maximum'].append(times[comp.maximum])
        wls = comp.created_from[0].dataset.wavelength
        spectrum_maxima, _ = find_peaks(comp.spectrum)
        spectrum_maxima = [m for m in spectrum_maxima if comp.spectrum[m] > 1]
        lambda_max = [wls[i] for i in spectrum_maxima]
        quali_comp_dict['lambda_max'].append(lambda_max)
        quali_comp_dict['num_peaks'].append(len(comp.created_from))
    return quali_comp_dict


def comps_to_df(comps):
    comp_dict = quali_comps_to_dict(comps)
    comp_df = pd.DataFrame(comp_dict)
    return comp_df


def create_quali_comp_page(comp):
    wls = comp.created_from[0].dataset.wavelength
    df = pd.DataFrame({'x': wls, 'y': comp.spectrum})
    spectrum = plot_1D_data(df, xlabel='Wavelength (nm)',
                            ylabel='Mean absorbance (mAU)',
                            title='',
                            reduce_data=False)

    peaks_df = peaks_to_df(comp.created_from)
    return dp.Page(
        title=comp.compound_id,
        blocks=[
            dp.Group(
                dp.Text(f"## Details to component {comp.compound_id}"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text(f"### Figure: UV-Vis spectrum of component {comp.compound_id}."),
            dp.Plot(spectrum),
            dp.Text("### Table: Peaks from which the component was created."),
            dp.DataTable(peaks_df, label=f"peak_table_{comp.compound_id}")
        ],        
    )


def report_quali_comps(quali_comp_db, report_path):
    comps = quali_comp_db.items
    comp_df = comps_to_df(comps)
    table_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# 3 Qualitative component database report"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text("### Table: Components in the qualitative component database."
                    "of the campaign."),
            dp.DataTable(comp_df, label="comp_table")
        ],        
    )
    component_pages = []
    for comp in quali_comp_db:
        page = create_quali_comp_page(comp)
        component_pages.append(page)
    r = dp.Report(
        table_page,
        *component_pages
    )
    r.save(path=os.path.join(report_path, "report_quali_comp_db.html"), open=True)
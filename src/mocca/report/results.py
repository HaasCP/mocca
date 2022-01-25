#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:08:02 2022

@author: haascp
"""

import os
import pandas as pd
import datapane as dp

from mocca.visualization.basic_plots import plot_1D_data


def chroms_to_results(chroms, quali_comp_db):
    compound_comps = [comp.compound_id for comp in quali_comp_db if not 
                      "unknown" in comp.compound_id and not "impurity" in comp.compound_id]
    
    chrom_dict = {'index': []}
    for comp in compound_comps:
        chrom_dict["integral_" + comp] = []
        chrom_dict["conc_" + comp] = []
    chrom_dict.update({'file': [],
                       'bad_data': [],
                       'compound_run': [],
                       'istd_added': []})

    
    for i, chrom in enumerate(chroms):
        chrom_dict['index'].append(i + 1)
        for comp in compound_comps:
            if any(peak.compound_id == comp for peak in chrom):
                chrom_dict["integral_" + comp].append(chrom[comp].integral)
                chrom_dict["conc_" + comp].append(chrom[comp].concentration)
            else:
                chrom_dict["integral_" + comp].append(None)
                chrom_dict["conc_" + comp].append(None)
        chrom_dict['file'].append(os.path.basename(chrom.dataset.path))
        chrom_dict['bad_data'].append(chrom.bad_data)
        chrom_dict['compound_run'].append(chrom.dataset.experiment.compound is not None)
        chrom_dict['istd_added'].append(chrom.dataset.experiment.istd is not None)
    return pd.DataFrame(chrom_dict)


def create_comp_pages(chroms, quali_comp_db, quant_comp_db):
    integral_dict = {comp.compound_id: [] for comp in quali_comp_db}
    conc_dict = {comp.compound_id: [] for comp in quant_comp_db}
    for chrom in chroms:
        for key in integral_dict.keys():
            if key in chrom:
                integral_dict[key].append(chrom[key].integral)
            else:
                integral_dict[key].append(0)
        for key in conc_dict.keys():
            if key in chrom:
                conc_dict[key].append(chrom[key].concentration)
            else:
                conc_dict[key].append(0)

    chrom_idxs = [i + 1 for i in list(range(len(chroms)))]
    comp_pages = []
    for key, val in integral_dict.items():
        df = pd.DataFrame({'idx': chrom_idxs,
                           key: val})
        comp_plot = plot_1D_data(df, xlabel='Chromatogram index', ylabel='Integral (mAU s)',
                                 title='', reduce_data=True)
        
        if key in conc_dict:
            df = pd.DataFrame({'idx': chrom_idxs,
                               key: conc_dict[key]})
            conc_plot = plot_1D_data(df, xlabel='Chromatogram index', ylabel='Concentration (mM)',
                                     title='', reduce_data=True)
        else:
            conc_plot = None
        
        blocks = [
            dp.Group(
                dp.Text(f"## Component {key} over runs"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text(f"### Figure: Integral of component {key} over runs."),
            dp.Plot(comp_plot)
        ]
        
        if conc_plot:
            blocks + [
                dp.Text(f"### Figure: Concentration of component {key} over runs."),
                dp.Plot(conc_plot)
                ]
        
        comp_page = dp.Page(
            title=key,
            blocks=blocks,
        )
        comp_pages.append(comp_page)
    return comp_pages


def report_runs(chroms, quali_comp_db, quant_comp_db, report_path):
    chrom_df = chroms_to_results(chroms, quali_comp_db)
    summary_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# 7 Results of compounds over runs"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text("### Table: Chromatograms processed during the campaign."),
            dp.DataTable(chrom_df, label="chrom_table")
        ],
    )
    comp_pages = create_comp_pages(chroms, quali_comp_db, quant_comp_db)
    r = dp.Report(
        summary_page,
        *comp_pages
    )
    r.save(path=os.path.join(report_path, "report_runs.html"), open=True)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:08:02 2022

@author: haascp
"""

import os
import pandas as pd
import datapane as dp

from mocca.visualization.basic_plots import plot_1D_scatter


def chroms_to_results(chroms, quali_comp_db):
    """
    Transfers relevant information from chromatograms in a pandas df.
    """
    compound_comps = [comp.compound_id for comp in quali_comp_db if
                      "unknown" not in comp.compound_id and
                      "impurity" not in comp.compound_id]

    chrom_dict = {'index': []}
    for comp in compound_comps:
        chrom_dict["integral_" + comp] = []
        chrom_dict["conc_" + comp] = []
    chrom_dict.update({'file': [],
                       'bad_data': [],
                       'compound_run': [],
                       'istd_added': []})

    for i, chrom in enumerate(chroms):
        if not chrom.bad_data:
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
            chrom_dict['compound_run'].append(chrom.experiment.compound is not None)
            chrom_dict['istd_added'].append(chrom.experiment.istd is not None)
    return pd.DataFrame(chrom_dict)


def create_comp_pages(chroms, quali_comp_db, quant_comp_db):
    """
    Creates a report page with details to each component over all runs of the
    cmapaign. Includes plots of integral and concentration over all runs.
    """
    integral_dict = {comp.compound_id: [] for comp in quali_comp_db}
    area_percent_dict = {comp.compound_id: [] for comp in quali_comp_db}
    conc_dict = {comp.compound_id: [] for comp in quant_comp_db}
    for chrom in chroms:
        if not chrom.bad_data:
            total_peak_sum = 0
            for peak in chrom.peaks:
                if peak.idx > 0:
                    total_peak_sum += peak.integral
            for key in integral_dict.keys():
                if key in chrom:
                    integral_dict[key].append(chrom[key].integral)
                    area_percent_dict[key].append(round(chrom[key].integral /
                                                        total_peak_sum * 100, 1))
                else:
                    integral_dict[key].append(0)
                    area_percent_dict[key].append(0)
            for key in conc_dict.keys():
                if key in chrom:
                    conc_dict[key].append(chrom[key].concentration)
                else:
                    conc_dict[key].append(0)

    chrom_idxs = [i + 1 for i in list(range(len(chroms))) if not chroms[i].bad_data]
    comp_pages = []
    for key, val in integral_dict.items():
        df = pd.DataFrame({'chromatogram_index': chrom_idxs,
                           key + '_integral': [round(i) for i in val]})
        comp_plot = plot_1D_scatter(df, xlabel='Chromatogram index',
                                    ylabel='Summed peak absorbance (mAU)',
                                    title='', reduce_data=True)

        df = pd.DataFrame({'chromatogram_index': chrom_idxs,
                           key + '_area_percent': area_percent_dict[key]})
        area_plot = plot_1D_scatter(df, xlabel='Chromatogram index',
                                    ylabel='Area percent (%)',
                                    title='', reduce_data=True)

        if key in conc_dict:
            df = pd.DataFrame({'chromatogram_index': chrom_idxs,
                               key + '_concentration': conc_dict[key]})
            conc_plot = plot_1D_scatter(df, xlabel='Chromatogram index',
                                        ylabel='Concentration (mM)',
                                        title='', reduce_data=True)
        else:
            conc_plot = None

        blocks = [
            dp.Group(
                dp.Text(f"## Component {key} over runs"),
                dp.Text("## MOCCA (Multivariate Online Contextual "
                        "Chromatographic Analysis)"),
                columns=2
            ),
            dp.Text(f"### Figure: Integral of component {key} over runs."),
            dp.Plot(comp_plot),
            dp.Text(f"### Figure: Area percent of component {key} over runs."),
            dp.Plot(area_plot)
        ]

        if conc_plot is not None:
            blocks = blocks + [
                dp.Text(f"### Figure: Concentration of component {key} over runs."),
                dp.Plot(conc_plot)
                ]

        comp_page = dp.Page(
            title=key,
            blocks=blocks,
        )
        comp_pages.append(comp_page)
    return comp_pages


def report_comp_tracking(chroms, quali_comp_db, quant_comp_db, report_path):
    """
    Main report function to follow concentrations and integrals over runs.
    """
    if not chroms or all([chrom.bad_data for chrom in chroms]):
        print("No chromatograms given or all chromatograms are bad data!")
        return
    this_dir, _ = os.path.split(__file__)
    mocca_icon_path = os.path.join(this_dir, "mocca_icon.png")
    chrom_df = chroms_to_results(chroms, quali_comp_db)
    summary_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# Results of compounds over runs"),
                dp.Media(file=mocca_icon_path),
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
    r.save(path=os.path.join(report_path, "compound_tracking.html"), open=True)

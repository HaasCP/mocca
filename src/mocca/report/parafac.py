#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 17:48:17 2022

@author: haascp
"""

import os
import math
import pandas as pd
import datapane as dp

from mocca.visualization.parafac_plots import plot_impure_peak_spectra, plot_parafac_peak_spec


def parafac_chroms_to_dict(chroms):
    chrom_dict = {'index': [],
                  'file': [],
                  'bad_data': [],
                  'compound_run': [],
                  'istd_added': [],
                  'num_peaks': []}

    for i, chrom in enumerate(chroms):
        if chrom.parafac_report_data:
            chrom_dict['index'].append(i + 1)
            chrom_dict['file'].append(os.path.basename(chrom.dataset.path))
            chrom_dict['bad_data'].append(chrom.bad_data)
            chrom_dict['compound_run'].append(chrom.dataset.experiment.compound is not None)
            chrom_dict['istd_added'].append(chrom.dataset.experiment.istd is not None)
            chrom_dict['num_peaks'].append(len(chrom.peaks))
    return chrom_dict


def parafac_chroms_to_df(chroms):
    chrom_dict = parafac_chroms_to_dict(chroms)
    chrom_df = pd.DataFrame(chrom_dict)
    return chrom_df


def create_parafac_page(chrom, index):
    
    for parafac_data in chrom.parafac_report_data:
        impure_peak = parafac_data[0]
        parafac_peaks = parafac_data[1]
        
        impure_peak_spec_plot = plot_impure_peak_spectra(impure_peak)
        parafac_spec_plots = [plot_parafac_peak_spec(peak) for peak in parafac_peaks]
        
        plots = [impure_peak_spec_plot] + parafac_spec_plots
        
        spectra = []
        for i in range(math.ceil(len(plots) / 2)):
            if 2 * i + 1 < len(plots):
                group = dp.Group(
                            dp.Plot(plots[2 * i]),
                            dp.Plot(plots[2 * i + 1]),
                            columns=2
                            )
            else:
                group = dp.Group(
                            dp.Plot(plots[2 * i]),
                            dp.Text("blank"),
                            columns=2
                            )
            spectra.append(group)

    return dp.Page(
        title=str(index),
        blocks=[
            dp.Group(
                dp.Text(f"## Details to chromatogram {index}"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text("### Figures: 1st figure is spectra of impure peak at all "
                    "time points. All other figures are spectra of parafac peaks "
                    "at elution maximum.")
        ] + spectra,      
    )


def report_parafac(chroms, report_path):
    chrom_df = parafac_chroms_to_df(chroms)
    summary_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# 7 PARAFAC report"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text("### Table: Chromatograms which triggered PARAFAC during the campaign."),
            dp.DataTable(chrom_df, label="chrom_table")
        ],
    )
    parafac_pages = []
    for i, chrom in enumerate(chroms):
        if chrom.parafac_report_data:
            page = create_parafac_page(chrom, i + 1)
            parafac_pages.append(page)
    r = dp.Report(
        summary_page,
        *parafac_pages
    )
    r.save(path=os.path.join(report_path, "report_parafac.html"), open=True)
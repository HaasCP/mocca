#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 17:48:17 2022

@author: haascp
"""

import os
import pandas as pd
import datapane as dp

from mocca.visualization.parafac_plots import (plot_impure_peak_spectra,
                                               plot_parafac_peaks_spec,
                                               plot_retention,
                                               plot_normalized_spectra,
                                               plot_normalized_elution,
                                               plot_normalized_integrals)


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
            chrom_dict['compound_run'].append(chrom.experiment.compound is not None)
            chrom_dict['istd_added'].append(chrom.experiment.istd is not None)
            chrom_dict['num_peaks'].append(len(chrom.peaks))
    return chrom_dict


def parafac_chroms_to_df(chroms):
    chrom_dict = parafac_chroms_to_dict(chroms)
    chrom_df = pd.DataFrame(chrom_dict)
    return chrom_df


def create_parafac_pages(chrom, index):
    
    parafac_pages = []
    for parafac_data in chrom.parafac_report_data:
        impure_peak = parafac_data[0]
        parafac_peaks = parafac_data[1]
        parafac_factors = parafac_data[2]
        
        impure_peak_spec_plot = plot_impure_peak_spectra(impure_peak)
        parafac_spec_plots = plot_parafac_peaks_spec(parafac_peaks)
        
        spectra = [impure_peak_spec_plot] + parafac_spec_plots
        
        retention_plot = plot_retention(impure_peak, parafac_peaks)
        
        normalized_spectra, normalized_elution, normalized_integrals = parafac_factors
        normalized_spectra_plot = plot_normalized_spectra(normalized_spectra)
        normalized_elution_plot = plot_normalized_elution(normalized_elution)
        normalized_integrals_plot = plot_normalized_integrals(normalized_integrals)
        

        parafac_page = dp.Page(
            title=f"chrom {str(index)}, peak {impure_peak.idx}",
            blocks=[
                dp.Group(
                    dp.Text(f"## Details to chromatogram {index}"),
                    dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                    columns=2
                ),
                dp.Text("### Figures: 1st figure is spectra of impure peak at all "
                        "time points. All other figures are spectra of parafac peaks "
                        "at elution maximum."),
                dp.Group(
                    *spectra,
                    columns=2
                    ),
                dp.Text("### Figure: Retention of impure peak overlayed with "
                        "calculated PARAFAC retention profiles."),
                dp.Plot(retention_plot),
                dp.Text("### Figure: PARAFAC model spectra of components."),
                dp.Plot(normalized_spectra_plot),
                dp.Text("### Figure: PARAFAC model elution profile of components."),
                dp.Plot(normalized_elution_plot),
                dp.Text("### Figure: RPARAFAC model integrals of components over runs."),
                dp.Plot(normalized_integrals_plot)
                ]
        )
        parafac_pages.append(parafac_page)
        return parafac_pages


def report_parafac(chroms, report_path):
    chrom_df = parafac_chroms_to_df(chroms)
    if chrom_df.empty:
        return
    summary_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# PARAFAC report"),
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
            pages = create_parafac_pages(chrom, i + 1)
            parafac_pages.extend(pages)
    r = dp.Report(
        summary_page,
        *parafac_pages
    )
    r.save(path=os.path.join(report_path, "report_parafac.html"), open=True)
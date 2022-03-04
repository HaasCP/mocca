#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 17:48:17 2022

@author: haascp
"""

import os
import pandas as pd
import datapane as dp

from mocca.visualization.parafac_plots import (plot_retention,
                                               plot_impure_peak_spectra,
                                               plot_uvvis_specs,
                                               plot_aligned_tensor,
                                               plot_normalized_integrals,
                                               plot_objective_func)


def parafac_chroms_to_df(chroms):
    """
    Transfers relevant information from chromatograms in a pandas df.
    """
    chrom_dict = {'index': [],
                  'file': [],
                  'bad_data': [],
                  'compound_run': [],
                  'istd_added': [],
                  'num_peaks': []}

    for i, chrom in enumerate(chroms):
        if chrom.parafac_models:
            chrom_dict['index'].append(i + 1)
            chrom_dict['file'].append(os.path.basename(chrom.dataset.path))
            chrom_dict['bad_data'].append(chrom.bad_data)
            chrom_dict['compound_run'].append(chrom.experiment.compound is not None)
            chrom_dict['istd_added'].append(chrom.experiment.istd is not None)
            chrom_dict['num_peaks'].append(len(chrom.peaks))
    chrom_df = pd.DataFrame(chrom_dict)
    return chrom_df


def create_parafac_pages(chrom, index):
    """
    Creates a report page with details to the PARAFAC model of a given impure
    peak. Only creates pages for models in which a known component was found.
    """
    parafac_pages = []
    report_models = [model for model in chrom.parafac_models if model.peaks]
    for parafac_model in report_models:
        retention_plot = plot_retention(parafac_model)

        impure_peak_spec_plot = plot_impure_peak_spectra(parafac_model.impure_peak)
        spec_plot = plot_uvvis_specs(parafac_model)

        aligned_retention_plot = plot_aligned_tensor(parafac_model)

        normalized_spectra, normalized_elution, normalized_integrals =\
            parafac_model.factors
        normalized_integrals_plot = plot_normalized_integrals(normalized_integrals)

        obj_func_plot = plot_objective_func(parafac_model)

        parafac_page = dp.Page(
            title=f"chrom {str(index)}, peak {parafac_model.impure_peak.idx}",
            blocks=[
                dp.Group(
                    dp.Text(f"## Details to chromatogram {index}, peak "
                            f"{parafac_model.impure_peak.idx}"),
                    dp.Text("## MOCCA (Multiway Online Chromatographic Chemical "
                            "Analysis)"),
                    columns=2
                ),
                dp.Text("### Figure: UV-Vis spectra of the impure peak at every "
                        "time point."),
                dp.Plot(impure_peak_spec_plot),
                dp.Text("### Figures: Modelled PARAFAC peaks after iterative "
                        "alignment."),
                dp.Group(
                    dp.Plot(retention_plot),
                    dp.Plot(spec_plot),
                    columns=2
                    ),
                dp.Text("### Figures: Visualization of the data tensor used for "
                        "the PARAFAC model."),
                dp.Group(
                    dp.Plot(aligned_retention_plot),
                    dp.Plot(normalized_integrals_plot),
                    columns=2
                    ),
                dp.Text("### Figure: Iterative PARAFAC objective function vs the "
                        "offset of the impure signal compared to the aligned pure "
                        "signals in the data tensor. Currently, the objective "
                        "function calculates the summed integral of all "
                        "components in the impure peak slice of the PARAFAC "
                        "model, which should be maximized."),
                dp.Plot(obj_func_plot)
                ]
        )
        parafac_pages.append(parafac_page)
        return parafac_pages


def report_parafac(chroms, report_path):
    """
    Main PARAFAC report function.
    """
    this_dir, _ = os.path.split(__file__)
    mocca_icon_path = os.path.join(this_dir, "mocca_icon.png")
    chrom_df = parafac_chroms_to_df(chroms)
    if chrom_df.empty:
        return
    summary_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# PARAFAC report"),
                dp.Media(file=mocca_icon_path),
                columns=2
            ),
            dp.Text("### Table: Chromatograms which triggered PARAFAC during "
                    "the campaign."),
            dp.DataTable(chrom_df, label="chrom_table")
        ],
    )
    parafac_pages = []
    for i, chrom in enumerate(chroms):
        if chrom.parafac_models:
            pages = create_parafac_pages(chrom, i + 1)
            if pages:
                parafac_pages.extend(pages)
    r = dp.Report(
        summary_page,
        *parafac_pages
    )
    r.save(path=os.path.join(report_path, "report_parafac.html"), open=True)

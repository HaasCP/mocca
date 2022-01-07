#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:12:43 2022

@author: haascp
"""
import os
import pandas as pd
import datapane as dp

from mocca.visualization.results_plot import plot_chrom_with_peaks
from mocca.report.utils import settings_to_df


def chroms_to_dict(chroms):
    chrom_dict = {'index': [],
                  'file': [],
                  'bad_data': [],
                  'compound_run': [],
                  'istd_added': [],
                  'num_peaks': []}

    for i, chrom in enumerate(chroms):
        chrom_dict['index'].append(i + 1)
        chrom_dict['file'].append(os.path.basename(chrom.dataset.path))
        chrom_dict['bad_data'].append(chrom.bad_data)
        chrom_dict['compound_run'].append(chrom.dataset.experiment.compound is not None)
        chrom_dict['istd_added'].append(chrom.dataset.experiment.istd is not None)
        chrom_dict['num_peaks'].append(len(chrom.peaks))
    return chrom_dict


def chroms_to_df(chroms):
    chrom_dict = chroms_to_dict(chroms)
    chrom_df = pd.DataFrame(chrom_dict)
    return chrom_df


def peaks_to_result_df(peaks):
    peaks_dict = {'peak_id': [],
                  'retention_time': [],
                  'compound_id': [],
                  'concentration': [],
                  'is_pure': [],
                  'is_saturated': [],
                  'integral': [],
                  'is_compound': []}
    for peak in peaks:
        times = peak.dataset.time
        peaks_dict['peak_id'].append(peak.idx)
        peaks_dict['retention_time'].append(times[peak.maximum + peak.offset])
        peaks_dict['compound_id'].append(peak.compound_id)
        peaks_dict['concentration'].append(peak.concentration)
        peaks_dict['is_pure'].append(peak.pure)
        peaks_dict['is_saturated'].append(peak.saturation)
        peaks_dict['integral'].append(peak.integral)
        peaks_dict['is_compound'].append(peak.is_compound)
    return pd.DataFrame(peaks_dict)

def create_chrom_page(chrom, index):
    chrom_plot = plot_chrom_with_peaks(chrom)

    peaks_df = peaks_to_result_df(chrom.peaks)
    return dp.Page(
        title=str(index),
        blocks=[
            dp.Group(
                dp.Text(f"## Details to chromatogram {index}"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text("### Figure: Chromatogram with highlighted peaks."),
            dp.Plot(chrom_plot),
            dp.Text("### Table: Peaks found in the chromatogram."),
            dp.DataTable(peaks_df)
        ],        
    )


def report_chroms(chroms, settings, report_path):
    chrom_df = chroms_to_df(chroms)
    summary_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# 6 Results by chromatogram report"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2
            ),
            dp.Text("### Table: Settings and thresholds used to process chromatograms."),
            dp.Table(settings_to_df(settings), label="settings_table"),
            dp.Text("### Table: Chromatograms processed during the campaign."),
            dp.DataTable(chrom_df, label="chrom_table")
        ],
    )
    chrom_pages = []
    for i, chrom in enumerate(chroms):
        page = create_chrom_page(chrom, i + 1)
        chrom_pages.append(page)
    r = dp.Report(
        summary_page,
        *chrom_pages
    )
    r.save(path=os.path.join(report_path, "report_chroms.html"), open=True)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 13:23:54 2022

@author: haascp
"""

import os
import pandas as pd
import datapane as dp

from mocca.report.hplc_input import exps_to_df
from mocca.dad_data.utils import sum_absorbance_by_time
from mocca.visualization.basic_plots import plot_1D_data


def report_gradient(gradient, report_path):
    exps = [gradient.experiment]
    exp_df = exps_to_df(exps)

    gradient_df = pd.DataFrame({
        'time': gradient.time,
        'absorbance': sum_absorbance_by_time(gradient.original_data)
        })
    gradient_plot = plot_1D_data(gradient_df, xlabel='Time (min)',
                                 ylabel='Summed absorbance (mAU s)',
                                 title='', reduce_data=True)
    baseline_df = pd.DataFrame({
        'time': gradient.time,
        'absorbance': sum_absorbance_by_time(gradient.data)
        })
    baseline_plot = plot_1D_data(baseline_df, xlabel='Time (min)',
                                 ylabel='Summed absorbance (mAU s)',
                                 title='', reduce_data=True)

    grad_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# Gradient report"),
                dp.Text("## MOCCA (Multiway Online Chromatographic Chemical Analysis)"),
                columns=2),
            dp.Text("### Table: Gradient experiment as given by the user."),
            dp.Table(exp_df, label="experiment_table"),
            dp.Group(
                dp.Text("### Plot: Original gradient."),
                dp.Text("### Plot: Baseline of gradient."),
                columns=2),
            dp.Group(
                dp.Plot(gradient_plot),
                dp.Plot(baseline_plot),
                columns=2)
        ],
    )
    r = dp.Report(
        grad_page
    )
    r.save(path=os.path.join(report_path, "report_gradient.html"), open=True)



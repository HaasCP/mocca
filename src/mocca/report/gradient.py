#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 13:23:54 2022

@author: haascp
"""

import os
import pandas as pd
import datapane as dp

from mocca.dad_data.utils import sum_absorbance_by_time
from mocca.visualization.basic_plots import plot_1D_data


def gradients_to_df(gradients):
    """
    Transfers relevant information of Gradient objects in a pandas df.
    """
    grad_dict = {'index': [],
                 'file': [],
                 'n_time_points': [],
                 'method_time': [],
                 'lambda_min': [],
                 'lambda_max': []}

    for i, grad in enumerate(gradients):
        grad_dict['index'].append(i + 1)
        grad_dict['file'].append(os.path.basename(grad.path))
        grad_dict['n_time_points'].append(len(grad.time))
        grad_dict['method_time'].append(round(max(grad.time), 2))
        grad_dict['lambda_min'].append(min(grad.wavelength))
        grad_dict['lambda_max'].append(max(grad.wavelength))
    grad_df = pd.DataFrame(grad_dict)
    return grad_df


def create_gradient_pages(gradients):
    """
    Creates a report page for each gradient used in the given campaign.
    """
    grad_pages = []
    for i, gradient in enumerate(gradients):
        gradient_df = pd.DataFrame({
            'time': gradient.time,
            'absorbance': sum_absorbance_by_time(gradient.original_data)
            })
        gradient_plot = plot_1D_data(gradient_df, xlabel='Time (min)',
                                     ylabel='Summed absorbance (mAU)',
                                     title='', reduce_data=True)
        baseline_df = pd.DataFrame({
            'time': gradient.time,
            'absorbance': sum_absorbance_by_time(gradient.data)
            })
        baseline_plot = plot_1D_data(baseline_df, xlabel='Time (min)',
                                     ylabel='Summed absorbance (mAU)',
                                     title='', reduce_data=True)
        exp_df = gradients_to_df([gradient])

        grad_page = dp.Page(
            title=f"Gradient {i + 1}",
            blocks=[
                dp.Text("### Table: Gradient details."),
                dp.Table(exp_df),
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
        grad_pages.append(grad_page)
    return grad_pages


def report_gradients(exps, report_path):
    """
    Main gradient report function.
    """
    if not any(exp.gradient for exp in exps):
        print("No gradients in the given HPLC input!")
        return
    this_dir, _ = os.path.split(__file__)
    mocca_icon_path = os.path.join(this_dir, "mocca_icon.png")
    gradients = []
    for exp in exps:
        if exp.gradient.dataset not in gradients:
            gradients.append(exp.gradient.dataset)
    grad_df = gradients_to_df(gradients)
    grad_pages = create_gradient_pages(gradients)

    summary_page = dp.Page(
        title="Start page",
        blocks=[
            dp.Group(
                dp.Text("# Gradient report"),
                dp.Media(file=mocca_icon_path),
                columns=2),
            dp.Text("### Table: Details to all gradients used in the campaign."),
            dp.Table(grad_df)
        ],
    )
    r = dp.Report(
        summary_page,
        *grad_pages
    )
    r.save(path=os.path.join(report_path, "report_gradient.html"), open=True)

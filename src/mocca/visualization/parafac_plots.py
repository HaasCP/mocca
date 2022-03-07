#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 18:05:10 2022

@author: haascp
"""
import altair as alt
import pandas as pd

from mocca.peak.utils import get_peak_data
from mocca.visualization.basic_plots import (plot_1D_layer, plot_1D_scatter_layer,
                                             plot_1D_scatter)
from mocca.dad_data.utils import sum_absorbance_by_time


color_vals = ['red', 'steelblue', 'chartreuse', '#F4D03F', '#D35400', '#7D3C98']


def plot_impure_peak_spectra(impure_peak):
    """
    Generates plot with UV-Vis spectra at every time point in the impure peak.
    """
    peak_data = get_peak_data(impure_peak)
    wls = impure_peak.dataset.wavelength

    charts = []
    for i in range(peak_data.shape[1]):
        df = pd.DataFrame({'x': wls,
                           'y': peak_data[:, i],
                           'z': impure_peak.dataset.time[impure_peak.left + i]})
        chart = alt.Chart(df, title='').mark_line().encode(
            x=alt.X(df.columns[0], axis=alt.Axis(title='Wavelength (nm)')),
            y=alt.Y(df.columns[1], axis=alt.Axis(title='Absorbance (mAU)')),
            color=alt.Color(df.columns[2], title='Time (min)',
                            scale=alt.Scale(scheme='viridis')),
            strokeWidth=alt.value(0.5))
        charts.append(chart)

    fig = charts[0]
    for i in range(len(charts) - 1):
        fig = fig + charts[i + 1]

    fig = fig.configure_axis(
        grid=False,
        titleFontSize=16,
        titleFontWeight='normal'
    ).configure_view(
        strokeWidth=0
    ).interactive()
    return fig


def plot_uvvis_specs(parafac_model):
    """
    Generates plot of the UV-Vis traces of the PARAFAC component together with
    the UV-Vis spectrum of the known compound used for the data tensor.
    """
    wls = parafac_model.impure_peak.dataset.wavelength
    fig = None
    for i, spectrum in enumerate(list(zip(*parafac_model.factors[0]))):
        df = pd.DataFrame({'x': wls,
                           'y': spectrum})
        p = plot_1D_layer(df, xlabel='Wavelength (nm)',
                          ylabel='Absorbance (mAU)',
                          title='', color=color_vals[i], reduce_data=True)
        if fig is None:
            fig = p
        else:
            fig += p

    tensor = parafac_model.data_tensor
    comp_spec = [val / (max(tensor.relevant_comp.spectrum) /
                        parafac_model.factors[0].max())
                 for val in tensor.relevant_comp.spectrum]

    df = pd.DataFrame({'x': wls,
                       'y': comp_spec})

    p = alt.Chart(df).mark_line(strokeDash=[10, 10]).encode(
                      x=alt.X(df.columns[0], axis=alt.Axis(title='Wavelength (nm)')),
                      y=alt.Y(df.columns[1], axis=alt.Axis(title='Absorbance (mAU)')),
                      color=alt.value('black')
                      ).properties(
                          title={"text": ['UV-Vis spectrum of component '
                                          f'{tensor.relevant_comp.compound_id} '
                                          '(black, dashed)',
                                          'and the computed PARAFAC peaks (colors)']})
    fig += p

    fig = fig.configure_axis(
        grid=False,
        titleFontSize=16,
        titleFontWeight='normal'
    ).configure_view(
        strokeWidth=0
    ).interactive()
    return fig


def plot_retention(parafac_model):
    """
    Plots normalized retention profiles of the PARAFAC components as well as of
    the impure peak.
    """
    impure_peak = parafac_model.impure_peak
    times = impure_peak.dataset.time

    peak_data = get_peak_data(impure_peak)
    summed_peak_data = sum_absorbance_by_time(peak_data)

    df = pd.DataFrame({'x': times[impure_peak.left:impure_peak.right + 1],
                       'y': summed_peak_data})

    impure_chart = alt.Chart(df).mark_line().encode(
                                 x=alt.X(df.columns[0],
                                         axis=alt.Axis(title='Time (min)')),
                                 y=alt.Y(df.columns[1],
                                         axis=alt.Axis(title='Absorbance (mAU)')),
                                 color=alt.value("black")
    ).properties(
        title={"text": ['Retention profiles of impure signal (black)',
                        'and of the computed PARAFAC peaks (colors)']})

    parafac_peaks = parafac_model.peaks
    left = parafac_peaks[0].left
    right = parafac_peaks[0].right + 1

    charts = []
    for i, peak in enumerate(parafac_peaks):
        peak_data = get_peak_data(peak)
        summed_peak_data = sum_absorbance_by_time(peak_data)
        df = pd.DataFrame({'x': times[left:right],
                           'y': summed_peak_data})
        chart = alt.Chart(df, title='').mark_line().encode(
            x=alt.X(df.columns[0], axis=alt.Axis(title='Time (min)')),
            y=alt.Y(df.columns[1], axis=alt.Axis(title='Absorbance (mAU)')),
            color=alt.value(color_vals[i])
        )
        charts.append(chart)

    fig = impure_chart
    for chart in charts:
        fig = fig + chart

    fig = fig.configure_axis(
        grid=False,
        titleFontSize=16,
        titleFontWeight='normal'
    ).configure_view(
        strokeWidth=0
    ).interactive()
    return fig


def plot_aligned_tensor(parafac_model):
    """
    Plots retention profiles (summed absorbance over all wavelengths) of all
    slices in the data tensor on the optimized iteration shift.
    """
    fig = None
    n_slices = parafac_model.data_tensor.tensor.shape[2]
    peak_tips = {'peak_max': []}
    for i in range(n_slices):
        if i == n_slices - 1:
            color = 'black'
        else:
            color = 'grey'
        summed_retention = parafac_model.data_tensor.tensor[:, :, i].sum(axis=0)
        peak_tips['peak_max'].append(summed_retention.argmax())
        df = pd.DataFrame({'x': list(range(len(summed_retention))),
                           'y': summed_retention})
        p = plot_1D_layer(df, xlabel='Time points (\u2013)',
                          ylabel='Absorbance (mAU)',
                          title='', color=color, reduce_data=True)
        if fig is None:
            fig = p
        else:
            fig += p

    rules = alt.Chart(pd.DataFrame(peak_tips)).mark_rule(strokeDash=[5, 5]).encode(
      x='peak_max')
    fig += rules

    fig = fig.properties(
        title={"text": ['Aligned retention profiles of impure signal (black)',
                        'and of the pure peaks (grey) composing the data tensor']}).\
        configure_axis(
            grid=False,
            titleFontSize=16,
            titleFontWeight='normal'
            ).configure_view(
                strokeWidth=0
                ).interactive()
    return fig


def plot_normalized_integrals(normalized_integrals):
    """
    Plots integrals of all PARAFAC components over all slices of the tensor.
    """
    fig = None
    for i, integral in enumerate(list(zip(*normalized_integrals))):
        df = pd.DataFrame({'x': list(range(len(integral))),
                           'y': integral})
        p = plot_1D_scatter_layer(df, xlabel='Data tensor slice (\u2013)',
                                  ylabel='Normalized integral (\u2013)',
                                  title='', color=color_vals[i], reduce_data=True)
        if fig is None:
            fig = p
        else:
            fig += p
    fig = fig.properties(
        title={"text": ['Integral values of the PARAFAC components over the '
                        'layers of the data tensor']}).\
        configure_axis(
            grid=False, titleFontSize=16, titleFontWeight='normal'
            ).configure_view(
                strokeWidth=0
                ).interactive()
    return fig


def plot_objective_func(parafac_model):
    """
    Plots objective function outcome vs iteration offset when using the iterative
    PARAFAC approach.
    """
    df = pd.DataFrame({'x': list(zip(*parafac_model.iter_objective_func))[0],
                       'y': list(zip(*parafac_model.iter_objective_func))[1]})
    fig = plot_1D_scatter(df, xlabel='Offset in time series points (\u2013)',
                          ylabel='Normalized integral (\u2013)',
                          title='', color='black', reduce_data=True)
    return fig

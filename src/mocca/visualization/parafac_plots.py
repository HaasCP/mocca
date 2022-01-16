#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 18:05:10 2022

@author: haascp
"""
import altair as alt
import pandas as pd
import numpy as np

from mocca.peak.utils import get_peak_data
from mocca.visualization.basic_plots import plot_1D_data
from mocca.dad_data.utils import sum_absorbance_by_time

def plot_impure_peak_spectra(impure_peak):
    peak_data = get_peak_data(impure_peak)
    wls = impure_peak.dataset.wavelength
    
    charts = []
    for i in range(peak_data.shape[1]):
        df = pd.DataFrame({'x': wls,
                           'y': peak_data[:, i]})
        chart = alt.Chart(df, title='').mark_line().encode(
            x=alt.X(df.columns[0], axis=alt.Axis(title='Wavelength (nm)')),
            y=alt.Y(df.columns[1], axis=alt.Axis(title='Absorbance (mAU)'))
        )
        charts.append(chart)
    
    fig = charts[0]
    for i in range(len(charts) - 1):
        fig = fig + charts[i + 1]
    
    fig = fig.configure_axis(
        grid=False,
        titleFontSize = 16,
        titleFontWeight='normal'
    ).configure_view(
        strokeWidth=0
    ).interactive()
    return fig


def plot_parafac_peak_spec(parafac_peak):
    wls = parafac_peak.dataset.wavelength
    
    peak_data = get_peak_data(parafac_peak)
    max_loc = np.argmax(np.sum(peak_data, axis=0))
    
    df = pd.DataFrame({'x': wls,
                       'y': peak_data[:, max_loc]})
    
    return plot_1D_data(df, xlabel='Wavelength (nm)', ylabel='Absorbance (mAU)',
                        title='', reduce_data=True)


def plot_retention(impure_peak, parafac_peaks):
    times = impure_peak.dataset.time
    peak_retention = sum_absorbance_by_time(impure_peak.dataset.data)
    
    left = parafac_peaks[0].left
    right = parafac_peaks[0].right
    
    df = pd.DataFrame({'x': times[(left + impure_peak.offset):\
                                  (right + impure_peak.offset + 1)],
                       'y': peak_retention[(left + impure_peak.offset):\
                                           (right + impure_peak.offset + 1)]})
    
    impure_chart = alt.Chart(df, title='').mark_line().encode(
        x=alt.X(df.columns[0], axis=alt.Axis(title='Time (min)')),
        y=alt.Y(df.columns[1], axis=alt.Axis(title='Absorbance (mAU)'))
    )
    
    charts = []
    for peak in parafac_peaks:
        summed_data = sum_absorbance_by_time(peak.dataset.data)
        df = pd.DataFrame({'x': times[left:right + 1],
                           'y': summed_data[left:right + 1]})
        chart = alt.Chart(df, title='').mark_line().encode(
            x=alt.X(df.columns[0], axis=alt.Axis(title='Time (min)')),
            y=alt.Y(df.columns[1], axis=alt.Axis(title='Absorbance (mAU)'))
        )
        charts.append(chart)
    
    fig = impure_chart
    for chart in charts:
        fig = fig + chart
    
    fig = fig.configure_axis(
        grid=False,
        titleFontSize = 16,
        titleFontWeight='normal'
    ).configure_view(
        strokeWidth=0
    ).interactive()
    return fig

    
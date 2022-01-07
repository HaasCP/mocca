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
    
    
    
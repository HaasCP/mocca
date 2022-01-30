#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 18:57:16 2021

@author: haascp
"""
import altair as alt
import numpy as np
import plotly.graph_objects as go


def plot_1D_data(df, xlabel='', ylabel='', title='', reduce_data=True):
    """
    Plots a set of 1D data.
    """
    if reduce_data:
        fac = df.shape[0] // 1000
        if fac > 0:
            df = df[::fac]

    chart = alt.Chart(df, title=title).mark_line().encode(
        x=alt.X(df.columns[0], axis=alt.Axis(title=xlabel)),
        y=alt.Y(df.columns[1], axis=alt.Axis(title=ylabel))
    ).configure_axis(
        grid=False,
        titleFontSize = 16,
        titleFontWeight='normal'
    ).configure_view(
        strokeWidth=0
    ).interactive()
    return chart


def plot_1D_scatter(df, xlabel='', ylabel='', title='', reduce_data=True):
    """
    Plots a set of 1D data.
    """
    if reduce_data:
        fac = df.shape[0] // 1000
        if fac > 0:
            df = df[::fac]

    chart = alt.Chart(df, title=title).mark_circle(size=60).encode(
        x=alt.X(df.columns[0], axis=alt.Axis(title=xlabel)),
        y=alt.Y(df.columns[1], axis=alt.Axis(title=ylabel)),
        tooltip=[df.columns[0], df.columns[1]]
    ).configure_axis(
        grid=False,
        titleFontSize = 16,
        titleFontWeight='normal'
    ).configure_view(
        strokeWidth=0
    ).interactive()
    return chart


def contour_map(data, time, wavelength, xlabel='Time',
                ylabel='Wavelength', title='Contour Plot', reduce_data=True):
    """
    Plots a contour map of the data.

    Parameters
    ----------
    data : numpy.array
        A numpy array of size (wavelength) x (timepoints) containing the
        data to be plotted

    time : list or numpy.array
        A list of timepoints to be plotted on the x-axis.

    wavelength : list or numpy.array
        A list of wavelengths to be plotted on the y-axis.

    xlabel : string, optional
        The x label of the output graph.
        Default is 'Time'

    ylabel : string, optional
        The y label of the output graph.
        Default is 'Wavelength'

    title : string, optional
        The title of the output graph.
        Default is 'Contour Plot'

    reduce_data : boolean
        If true, then the number of data points in the visualization is lowered.
        Default is True.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        A contour plot of the data.

    """
    if reduce_data:
        wl_fac = 1  # max(1, data.shape[0] // 100)
        t_fac = max(1, data.shape[1] // 1000)
        time = time[::t_fac]
        wavelength = wavelength[::wl_fac]
        data = data[::wl_fac, ::t_fac]

    # normalize absorbance and remove values close to 0
    data = np.maximum(data, 0) / data.max()
    data[data < 0.01] = 0

    fig = go.Figure(data=go.Contour(
            z=data,
            x=time,
            y=wavelength,
            colorbar=dict(
                title='Normalized absorbance (\u2013)',  # title here
                titleside='right',
                nticks=10
            ),
            contours_coloring='lines',
            colorscale='viridis',
            contours=dict(
                start=0,
                end=1,
                size=0.05,
            ),
            line_width=1
        ))
    fig.update_layout(
        title=title,
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        autosize=False,
        width=750,
        height=500)
    return fig

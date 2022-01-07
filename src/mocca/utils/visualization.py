#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 18:57:16 2021

@author: haascp
"""
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import altair as alt


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
        wl_fac = max(1, data.shape[0] // 100)
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


def plot_1D_data(df, xlabel='', ylabel='', title='', reduce_data=True):
    """
    Plots a set of 1D data at a specific wavelength against time.

    Parameters
    ----------
    df : pandas.DataFrame
        df should be a 2-column dataframe, consisting of the x-variable
        in the 0th column and the y-variable in the 1st column.

    wavelength : list or numpy.array
        A list of wavelengths to be plotted on the y-axis.

    reduce_data : boolean
        If true, then the number of data points in the visualization is lowered.
        Default is True.

    xlabel : string, optional
        The x label of the output graph.
        Default is ''.

    ylabel : string, optional
        The y label of the output graph.
        Default is ''.

    title : string, optional
        The title of the output graph.
        Default is ''.

    reduce_data : boolean
        If true, then the number of data points in the visualization is lowered.
        Default is True.

    Returns
    -------
    chart : altair.Chart
        A plot of the data.

    """
    if reduce_data:
        fac = df.shape[0] // 1000
        if fac > 0:
            df = df[::fac]

    # max_absorbance = df[df.columns[1]].max()
    # df.loc[:, df.columns[1]] = df[df.columns[1]] / max_absorbance

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


def round_to_n(x, n):
    """
    Returns number in a format suitable for data visualization.
    """
    if 1e-3 < x < 1e3:
        return round(x, n)
    else:
        return "{0:.{1}e}".format(x, n)


def plot_calibration_curves(comp):
    calibration_curves = []
    for version in comp.calib_factors.keys():
        # visualization
        curve_annotation_formula = f"y = {round_to_n(comp.calib_factors[version], 2)} x"
        curve_annotation_accuracy = f"R\u00B2 = {round(comp.calib_scores[version], 4)}"
        title = f"{version}"

        xlabel='Concentration (M)'
        if version != 'absolute':           
            ylabel=f"Relative area * {version} concentration (M)"
        else:
            ylabel='Area (mAU s)'
        
        x, y = zip(*comp.calib_data[version])
        
        df_scatter = pd.DataFrame({
            'x': x,
            'y': y
        })
        df_line = pd.DataFrame({
            'x': x,
            'y': [item * comp.calib_factors[version] for item in x]
        })
    
        scatter = alt.Chart(df_scatter, title=title).mark_circle(size=80, opacity=1).encode(
            x=alt.X(df_scatter.columns[0], axis=alt.Axis(title=xlabel)),
            y=alt.Y(df_scatter.columns[1], axis=alt.Axis(title=ylabel)),
            tooltip=[alt.Tooltip(df_scatter.columns[0], title=xlabel), 
                     alt.Tooltip(df_scatter.columns[1], title=ylabel)]
        ).interactive()
        
        chart = alt.Chart(df_line).mark_line(color='black').encode(
                x=alt.X(df_line.columns[0], scale = alt.Scale(domain = (0.9 * np.min(x), np.max(x) * 1.1))),
                y=alt.Y(df_line.columns[1], scale = alt.Scale(domain = (0.9 * min(np.min(y), np.min([item / comp.calib_factors[version] for item in x])), 
                                                                        1.1 * max(np.max(y), np.max([item / comp.calib_factors[version] for item in x])))))
            )

        if len(df_scatter[df_scatter.columns[0]]) > 1:
            annotation_x_loc = (min(df_scatter[df_scatter.columns[0]]) +
                                (max(df_scatter[df_scatter.columns[0]] -
                                     min(df_scatter[df_scatter.columns[0]]))
                                 * 0.1))
            annotation_y_loc_1 = (min(df_scatter[df_scatter.columns[1]]) +
                                  (max(df_scatter[df_scatter.columns[1]] -
                                       min(df_scatter[df_scatter.columns[1]]))
                                   * 0.95))
            annotation_y_loc_2 = (min(df_scatter[df_scatter.columns[1]]) +
                                  (max(df_scatter[df_scatter.columns[1]] -
                                       min(df_scatter[df_scatter.columns[1]]))
                                   * 0.75))
        else:
            annotation_x_loc = df_scatter[df_scatter.columns[0]][0] * 0.98
            annotation_y_loc_1 = df_scatter[df_scatter.columns[1]][0] * 0.95
            annotation_y_loc_2 = df_scatter[df_scatter.columns[1]][0] * 0.9

        annotation_formula = alt.Chart({'values':[{'x': annotation_x_loc,
                                                   'y': annotation_y_loc_1}]}).mark_text(
            text=curve_annotation_formula, align='left'
        ).encode(
            x='x:Q', y='y:Q'
        )

        annotation_accuracy = alt.Chart({'values':[{'x': annotation_x_loc,
                                                   'y': annotation_y_loc_2}]}).mark_text(
            text=curve_annotation_accuracy, align='left'
        ).encode(
            x='x:Q', y='y:Q'
        )

        fig = chart + scatter + annotation_formula + annotation_accuracy
        fig = fig.configure_axis(
                grid=False,
                titleFontSize = 16,
                titleFontWeight='normal'
            ).configure_view(
                strokeWidth=0
            )
        calibration_curves.append(fig)
    return calibration_curves

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 18:57:16 2021

@author: haascp
"""
import altair as alt


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


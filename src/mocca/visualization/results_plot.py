#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:12:02 2022

@author: haascp
"""
import pandas as pd
import altair as alt

from mocca.dad_data.utils import sum_absorbance_by_time


def plot_chrom_with_peaks(chrom):
    """Plots summed absorbance vs time with highlighted picked peak zones."""

    df = pd.DataFrame({
        'time': chrom.dataset.time,
        'absorbance': sum_absorbance_by_time(chrom.dataset.data)
        })
    fac = df.shape[0] // 1000
    if fac > 0:
        df = df[::fac]

    xlabel = 'Time (min)'
    ylabel = 'Summed absorbance (mAU)'
    title = ''

    chart = alt.Chart(df, title=title).mark_line().encode(
        x=alt.X(df.columns[0], axis=alt.Axis(title=xlabel)),
        y=alt.Y(df.columns[1], axis=alt.Axis(title=ylabel))
    )

    peaks = [peak for peak in chrom if peak.idx > 0]
    areas = []
    borders = []
    for peak in peaks:
        # choose color
        if peak.pure and not peak.saturation:
            color = 'green'
        elif peak.pure and peak.saturation:
            color = 'orange'
        else:
            color = 'red'
        peak_df = pd.DataFrame({
            'start': [chrom.dataset.time[peak.left]],
            'stop': [chrom.dataset.time[peak.right]]
        })
        area = alt.Chart(
            peak_df
        ).mark_rect(
            opacity=0.2,
            color=color
        ).encode(
            x='start',
            x2='stop',
            y=alt.value(0),  # pixels from top
            y2=alt.value(300)  # pixels from top
        )
        areas.append(area)

        brdrs = pd.DataFrame({'peak_borders': [chrom.dataset.time[peak.left],
                                               chrom.dataset.time[peak.right]]})
        border = alt.Chart(brdrs).mark_rule(color=color).encode(
          x='peak_borders')
        borders.append(border)

    peak_tips = {'peak_max': []}
    for peak in chrom:
        peak_tips['peak_max'].append(chrom.dataset.time[peak.maximum])
    rules = alt.Chart(pd.DataFrame(peak_tips)).mark_rule(strokeDash=[5, 5]).encode(
      x='peak_max')

    fig = chart + rules
    for area in areas:
        fig = fig + area
    for border in borders:
        fig = fig + border

    fig = fig.configure_axis(
        grid=False,
        titleFontSize=16,
        titleFontWeight='normal'
    ).configure_view(
        strokeWidth=0
    ).interactive()
    return fig

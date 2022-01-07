#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 14:15:20 2022

@author: haascp
"""
import numpy as np
import pandas as pd
import altair as alt

from mocca.visualization.utils import round_to_n


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
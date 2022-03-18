#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 15:49:54 2022

@author: haascp
"""

from mocca.peak.models import ProcessedPeak


def quantify_peak(peak, quant_comp_db):
    """
    Takes a peak with an integral and quantifies it by calculating the corresponding
    concentration if the compound is present in the quantification database.
    """
    if peak.compound_id in quant_comp_db:
        quant_comp = quant_comp_db[peak.compound_id]
        if peak.istd and not any(istd_peak.compound_id == peak.compound_id for
                                 istd_peak in peak.istd):
            scores = quant_comp.calib_scores
            max_score_version = ''
            max_score = 0
            for istd_p in [peak for peak in peak.istd if peak.concentration]:
                if scores[istd_p.compound_id] > max_score:
                    max_score_version = istd_p.compound_id
                    max_score = scores[istd_p.compound_id]
                    istd_peak = istd_p
            concentration = (peak.integral * istd_peak.concentration /
                             istd_peak.integral /
                             quant_comp.calib_factors[max_score_version])

        else:
            concentration = (peak.integral /
                             quant_comp.calib_factors['absolute'])
    else:
        concentration = None

    return ProcessedPeak(left=peak.left,
                         right=peak.right,
                         maximum=peak.maximum,
                         dataset=peak.dataset,
                         idx=peak.idx,
                         saturation=peak.saturation,
                         pure=peak.pure,
                         integral=peak.integral,
                         istd=peak.istd,
                         offset=peak.offset,
                         compound_id=peak.compound_id,
                         concentration=concentration,
                         is_compound=peak.is_compound)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 15:49:54 2022

@author: haascp
"""

from mocca.peak.models import ProcessedPeak

def quantify_peak(peak, quant_comp_db):
    if peak.compound_id in quant_comp_db:
        quant_comp = quant_comp_db[peak.compound_id]
        if peak.istd:
            scores = quant_comp.calib_scores
            max_score_version = 'absolute'
            max_score = scores['absolute']
            for istd_p in [peak for peak in peak.istd if peak.concentration]:
                if scores[istd_p.compound_id] > max_score:
                    max_score_version = istd_p.compound_id
                    max_score = scores[istd_p.compound_id]
                    istd_peak = istd_p
            if max_score_version == 'absolute':
                concentration = (peak.integral /
                                 quant_comp.calib_factors['absolute'])
            else:
                concentration = (peak.integral * istd_peak.concentration / istd_peak.integral /
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
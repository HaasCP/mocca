#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 08:59:36 2022

@author: haascp
"""
import numpy as np

from mocca.peak.models import CorrectedPeak, IntegratedPeak
from mocca.dad_data.models import ParafacData


def create_pure_peak(impure_peak):
    """
    Takes an impure peak and returns its copy with the pure attribute True.
    """
    if type(impure_peak) == CorrectedPeak:
        parafac_peak = CorrectedPeak(left=impure_peak.left,
                                     right=impure_peak.right,
                                     maximum=impure_peak.maximum,
                                     offset=impure_peak.offset,
                                     dataset=impure_peak.dataset,
                                     idx=-impure_peak.idx,
                                     saturation=impure_peak.saturation,
                                     pure=True,
                                     integral=impure_peak.integral,
                                     istd=impure_peak.istd)
    elif type(impure_peak) == IntegratedPeak:
        parafac_peak = IntegratedPeak(left=impure_peak.left,
                                      right=impure_peak.right,
                                      maximum=impure_peak.maximum,
                                      offset=impure_peak.offset,
                                      dataset=impure_peak.dataset,
                                      idx=-impure_peak.idx,
                                      saturation=impure_peak.saturation,
                                      pure=True,
                                      integral=impure_peak.integral)
    else:
        raise TypeError(f"Given impure peak is of type {type(impure_peak)}. "
                        "Only mocca IntegratedPeak and CorrectedPeak types "
                        "are allowed.")
    return parafac_peak


def get_parafac_data_shift(iter_offset):
    """
    If the iteration offset is larger than zero, the impure peak was shifted
    and must therefore be shifted back in the resulting PARAFAC data. If the
    iteration shift was negative, the pure signals were shifted and nothing
    has to be done.
    """
    if iter_offset > 0:
        shift = iter_offset
    else:
        shift = 0
    return shift


def create_parafac_peak(comp_i, parafac_model):
    """
    Return synthetic PARAFAC peaks created from the PARAFAC decomposition results.
    """
    shift = get_parafac_data_shift(parafac_model.iter_offset)
    parafac_comp_factors = (parafac_model.factors[0][:, comp_i],
                            parafac_model.factors[1][:, comp_i],
                            parafac_model.factors[2][:, comp_i])
    tensor = parafac_model.data_tensor
    impure_peak = parafac_model.impure_peak
    left_bound = tensor.boundaries[0]
    right_bound = tensor.boundaries[1]

    if type(impure_peak) == CorrectedPeak:
        parafac_peak = CorrectedPeak(left=left_bound,
                                     right=right_bound,
                                     maximum=(left_bound - shift +
                                              np.argmax(parafac_comp_factors[1])),
                                     offset=impure_peak.offset,
                                     dataset=ParafacData(impure_peak,
                                                         parafac_comp_factors,
                                                         tensor.boundaries,
                                                         shift,
                                                         tensor.y_offset),
                                     idx=-impure_peak.idx,
                                     saturation=impure_peak.saturation,
                                     pure=True,
                                     # reaction run in last slice of data tensor
                                     integral=(parafac_comp_factors[2][-1]),
                                     istd=impure_peak.istd)
    elif type(impure_peak) == IntegratedPeak:
        parafac_peak = IntegratedPeak(left=left_bound,
                                      right=right_bound,
                                      maximum=(left_bound - shift +
                                               np.argmax(parafac_comp_factors[1]) -
                                               shift),
                                      offset=impure_peak.offset,
                                      dataset=ParafacData(impure_peak,
                                                          parafac_comp_factors,
                                                          tensor.boundaries,
                                                          shift,
                                                          tensor.y_offset),
                                      idx=-impure_peak.idx,
                                      saturation=impure_peak.saturation,
                                      pure=True,
                                      # reaction run in last slice of data tensor
                                      integral=(parafac_comp_factors[2][-1]))
    else:
        raise TypeError(f"Given impure peak is of type {type(impure_peak)}. "
                        "Only mocca IntegratedPeak and CorrectedPeak types "
                        "are allowed.")
    return parafac_peak

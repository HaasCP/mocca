#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:05:19 2021

@author: haascp
"""
import numpy as np

from mocca.peak_utils import get_peak_data
def integrate_peak(assigned_peak):
    """
    Integrates the peak.

    Modifies
    --------
    picked_peak.integral : Integral of the peak
    """
    peak_data = get_peak_data(assigned_peak)
    return np.sum(peak_data).tolist()

def quantify_peak(integrated_peak, quant_component_db):
    """
    Computes the concentration of the compound in this peak.

    The attributes self.integral and self.compound_id must be set beforehand
    (through functions self.integrate_peak() and self.check_database())
    in order to quantify.

    Parameters
    ----------
    quantification_database : micdrop.QuantificationDatabase
        A quantification database containing all seen components.

    Raises Exception if the attributes self.integral or self.compound_id are
    not set. Prints a text warning if self.pure is not set.

    Modifies
    --------
    self.concentration : sets concentration to that predicted by integral
    """
    if integrated_peak.compound_id is None:
        return None

    if integrated_peak.pure is True:
        pass
    else:
        return None

    self.concentration = quantification_database.quantify_peak(self.integral,
                                                               self.compound_id)
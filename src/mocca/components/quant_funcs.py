#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 17:20:35 2021

@author: haascp
"""
import numpy as np
from sklearn.linear_model import LinearRegression

from mocca.components.models import QuantComponent

def create_quant_component(peaks):
    """
    Creates a quantitative component object based on the given peaks
    """
    warnings = []
    if not peaks:
        return
    if all(peak.compound_id == peaks[0].compound_id for peak in peaks):
        compound_id = peaks[0].compound_id
    else:
        raise AttributeError("All peaks have to have the same compound_id to "
                             "create a component")
    calib_data = {}
    calib_data['absolute'] = []
    for peak in peaks:
        calib_point = (peak.concentration, peak.integral)
        calib_data['absolute'].append(calib_point)

        if peak.istd:
            for istd_peak in peak.istd:
                if not istd_peak.compound_id in calib_data:
                    calib_data[istd_peak.compound_id] = []
                calib_point = (peak.concentration, (peak.integral *
                                                    istd_peak.concentration /
                                                    istd_peak.integral))
                calib_data[istd_peak.compound_id].append(calib_point)
    
    calib_factors = {}
    for version, calibration_data in calib_data.items(): 
        calibration_curve = LinearRegression(fit_intercept = False)
        x, y = zip(*calibration_data)
        x = np.array(x).reshape((-1,1))
        calibration_curve.fit(x, y)
        score = calibration_curve.score(x, y)
        if score < 0.99:
            warnings.append("Calibration curve of version '{}' has low "
                            "correlation of {}!".format(version, score))
        
        calib_factors[version] = calibration_curve.coef_[0]

    return QuantComponent(compound_id=compound_id,
                          calib_factors=calib_factors,
                          calib_data=calib_data,
                          created_from=peaks)

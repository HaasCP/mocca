#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 17:20:35 2021

@author: haascp
"""
import numpy as np
from sklearn.linear_model import LinearRegression

from mocca.peak.utils import get_peak_data
from mocca.components.models import QuantComponent
from mocca.components.utils import check_peaks_compound_id


def get_integrate_wl_index(compound_id, quali_comp_db):
    """
    Returns index of the wavelength vector where spectrum of component has the
    highest maximum.
    """
    if quali_comp_db[compound_id].spectrum_max:
        integrate_wl_idx = quali_comp_db[compound_id].spectrum_max[0]
    else:
        integrate_wl_idx = 1  # minimum wavelength due to bandwidth of 2
    return integrate_wl_idx


def integrate_on_wl(peak, integrate_wl_idx, bandwidth=2):
    """
    Integrates signal of given peak on a given wavelength with a bandwidth of 2
    (default) by summing all absorbances.
    """
    peak_data = get_peak_data(peak)
    peak_data_wl = peak_data[integrate_wl_idx - 1:integrate_wl_idx + 2]
    # correct baseline
    peak_data_wl = peak_data_wl - peak_data_wl.min()

    integral = np.sum(peak_data_wl).tolist()
    return integral


def create_calibration_dict(peaks, integrate_wl_idx, quali_comp_db):
    """
    Creates a dictionary with all data needed to create calibration curves.
    """
    calib_data = {}
    calib_data['absolute'] = []
    for peak in peaks:
        integral_wl = integrate_on_wl(peak, integrate_wl_idx)
        calib_point = (peak.concentration, integral_wl)
        calib_data['absolute'].append(calib_point)
        if peak.istd:
            for istd_peak in peak.istd:
                istd_int_wl_index = get_integrate_wl_index(istd_peak.compound_id,
                                                           quali_comp_db)
                istd_int_wl = integrate_on_wl(istd_peak, istd_int_wl_index)

                if istd_peak.compound_id not in calib_data:
                    calib_data[istd_peak.compound_id] = []
                calib_point = (peak.concentration, (integral_wl *
                                                    istd_peak.concentration /
                                                    istd_int_wl))
                calib_data[istd_peak.compound_id].append(calib_point)
    return calib_data


def create_linear_models(calib_data):
    """
    Creates linear models out of the given data and returns a list of linear
    calibration factors (calibration curve is forced through origin) and a list
    of corresponing R-squared values of the regression.
    """
    calib_factors = {}
    calib_scores = {}
    for version, calibration_data in calib_data.items():
        calibration_curve = LinearRegression(fit_intercept=False)
        x, y = zip(*calibration_data)
        x = np.array(x).reshape((-1, 1))
        calibration_curve.fit(x, y)
        score = calibration_curve.score(x, y)

        calib_factors[version] = calibration_curve.coef_[0]
        calib_scores[version] = score
    return calib_factors, calib_scores


def create_quant_component(peaks, quali_comp_db):
    """
    Creates a quantitative component object based on the given peaks
    """
    if not peaks:
        return None

    compound_id = check_peaks_compound_id(peaks)

    integrate_wl_idx = get_integrate_wl_index(compound_id, quali_comp_db)

    calib_data = create_calibration_dict(peaks, integrate_wl_idx, quali_comp_db)

    calib_factors, calib_scores = create_linear_models(calib_data)

    return QuantComponent(compound_id=compound_id,
                          integrate_wl_idx=integrate_wl_idx,
                          calib_factors=calib_factors,
                          calib_data=calib_data,
                          calib_scores=calib_scores,
                          created_from=peaks)

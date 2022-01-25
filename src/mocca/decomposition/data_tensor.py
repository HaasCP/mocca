#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 09:04:46 2022

@author: haascp
"""
import numpy as np

from mocca.peak.utils import get_peak_data
from mocca.decomposition.utils import check_comp_overlap


def get_relevant_comps(impure_peak, quali_comp_db):
    """
    Returns components which are not unknown or impurities (-> only components
    which were given by the user via compound) and which overlap with the impure
    peak.
    """
    return [comp for comp in quali_comp_db if (check_comp_overlap(impure_peak, comp)
                                               and (not 'unknown' in comp.compound_id
                                                    and not 'impurity' in comp.compound_id))]


def get_offset_peak_to_comp(created_from_peak, comp):
    """
    Returns the time offset of the maxima of oeak and quali component
    """
    return created_from_peak.maximum - comp.maximum


def get_comp_peaks(relevant_comp):
    """
    Returns up to five peaks from which the given component was created.
    """
    created_from_peaks = relevant_comp.created_from
    # take maximum 5 peaks
    fac = len(created_from_peaks) // 5
    if fac > 0:
        created_from_peaks = created_from_peaks[::fac]
    return created_from_peaks


def get_tensor_boundaries(impure_peak, relevant_comps, iter_offset):
    """
    Returns a tuple containing the leftest and rightest peak/component boundaries.
    The boundaries are corrected by the iteration offset if the iterative
    PARAFAC algorithm is used.
    """
    lefts = []
    rights = []
    for comp in relevant_comps:
        created_from_peaks = get_comp_peaks(comp)
        for peak in created_from_peaks:
            offset = get_offset_peak_to_comp(peak, comp)
            lefts.append(peak.left - offset)
            rights.append(peak.right - offset)

    left_boundaries = lefts + [impure_peak.left]
    right_boundaries = rights + [impure_peak.right]

    left = min(left_boundaries)
    right = max(right_boundaries)

    if iter_offset < 0:
        right = right - iter_offset
    elif iter_offset > 0:
        right = right + iter_offset
    return (left, right)


def get_zeros_array(boundaries, n_wavelengths):
    """
    Retruns array of zeros in the shape of the wavelengths and boundaries.
    """
    return np.zeros((n_wavelengths, boundaries[1] - boundaries[0] + 1))


def get_zero_extended_peak_data(peak_data, left, boundaries):
    """
    Returns zero-extended (to the boundaries) peak data.
    """
    peak_data_ze = get_zeros_array(boundaries, peak_data.shape[0])
    rel_left = left - boundaries[0]
    rel_right = rel_left + peak_data.shape[1]
    peak_data_ze[:, rel_left:rel_right] = peak_data
    return peak_data_ze


def get_comp_peak_data_list(relevant_comps, boundaries, iter_offset):
    """
    Returns a list of maximum-aligned peak data of the comp peaks of the
    relevant components. Moreover, returns the shape of the component list
    (how many slices per component in the list).
    """
    ze_peaks = []
    comp_tensor_shape = tuple()
    for comp in relevant_comps:
        created_from_peaks = get_comp_peaks(comp)
        comp_tensor_shape = comp_tensor_shape + (len(created_from_peaks),)
        for peak in created_from_peaks:
            offset = get_offset_peak_to_comp(peak, comp)
            peak_data = get_peak_data(peak)
            if iter_offset < 0:
                left = peak.left - offset - iter_offset
            else:
                left = peak.left - offset
            peak_data_ze = get_zero_extended_peak_data(peak_data,
                                                       left,
                                                       boundaries)
            if peak_data_ze.min() < 0:
                peak_data_ze = peak_data_ze - peak_data_ze.min()
            ze_peaks.append(peak_data_ze)
    return ze_peaks, comp_tensor_shape


def get_zero_extended_impure_peak_data(impure_peak, boundaries, iter_offset):
    """
    Returns zero-extended (to the boundaries) peak data of the impure peak.
    The location of the peak is corrected by the iteration offset in case
    the iterative PARAFAC algorithm is used.
    """
    peak_data = get_peak_data(impure_peak)
    if iter_offset > 0:
        left = impure_peak.left + iter_offset
    else:
        left = impure_peak.left
    peak_data_ze = get_zero_extended_peak_data(peak_data,
                                               left,
                                               boundaries)
    if peak_data_ze.min() < 0:
        y_offset = peak_data_ze.min()
        peak_data_ze = peak_data_ze - peak_data_ze.min()
    else:
        y_offset = 0
    return peak_data_ze, y_offset

def get_parafac_data_list(impure_peak, quali_comp_db, iter_offset):
    """
    Returns a list of peak data which will be transferred into the PARAFAC
    data tensor
    """
    relevant_comps = get_relevant_comps(impure_peak, quali_comp_db)

    boundaries = get_tensor_boundaries(impure_peak, relevant_comps, iter_offset)

    comp_peaks, comp_tensor_shape = get_comp_peak_data_list(relevant_comps,
                                                            boundaries,
                                                            iter_offset)

    impure_peak, y_offset = get_zero_extended_impure_peak_data(impure_peak,
                                                               boundaries,
                                                               iter_offset)

    return (comp_peaks + [impure_peak], boundaries, relevant_comps,
            comp_tensor_shape, y_offset)


def get_parafac_tensor(impure_peak, quali_comp_db, iter_offset,
                       show_parafac_analytics):
    """
    Returns the data tensor used for subsequent PARAFAC decomposition.
    """
    parafac_data, boundaries, relevant_comps, comp_tensor_shape, y_offset = \
        get_parafac_data_list(impure_peak, quali_comp_db, iter_offset)

    # add all data to data_flattened and create data tensor
    data_flattened = []
    for data in parafac_data:
        data_flattened.append(data[:, :, np.newaxis])
    data_tensor = np.concatenate(data_flattened, axis=2)
    if show_parafac_analytics:
        print(f"new data tensor with boundaries {boundaries} and shape"
              f"{comp_tensor_shape} built from "
              f"{[comp.compound_id for comp in relevant_comps]}.")
    return data_tensor, boundaries, relevant_comps, comp_tensor_shape, y_offset

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 09:04:46 2022

@author: haascp
"""
import numpy as np

from mocca.dad_data.utils import sum_absorbance_by_time
from mocca.peak.utils import get_peak_data
from mocca.decomposition.utils import check_comp_overlap
from mocca.decomposition.model import DataTensor


def get_relevant_comp(impure_peak, quali_comp_db):
    """
    Returns component which is not unknown or impurity (-> only components
    which were given by the user via compound) and which overlap with the impure
    peak. If multiple components overlapping the component with the best
    UV-Vis correlation coefficient is returned.
    """
    relevant_comps = [comp for comp in quali_comp_db if
                      (check_comp_overlap(impure_peak, comp) and
                       ('unknown' not in comp.compound_id and
                        'impurity' not in comp.compound_id))]
    relevant_comp = sorted(relevant_comps,
                           key=lambda comp: np.corrcoef(
                               impure_peak.dataset.data[:, impure_peak.maximum],
                               comp.spectrum)[0, 1])[-1]
    return relevant_comp


def get_offset_peak_to_comp(created_from_peak, comp):
    """
    Returns the time offset of the maxima of oeak and quali component.
    """
    return created_from_peak.maximum - comp.maximum


def get_comp_peaks(relevant_comp):
    """
    Returns exactly five peaks from which the given component was created. If
    less than 5 peaks were used to create the component, the peaks are
    multiplied as long as five peaks are reached.
    """
    created_from_peaks = relevant_comp.created_from
    # take maximum 2 peaks
    fac = len(created_from_peaks) // 2
    if fac > 2:
        created_from_peaks = created_from_peaks[::fac]
    # duplicate peaks to have at least 2 peaks
    elif len(created_from_peaks) < 2:
        num_repeats = int(2 / len(created_from_peaks))
        new_created_from_peaks = [val for val in created_from_peaks for
                                  _ in range(num_repeats)]
        i = 0
        while len(new_created_from_peaks) < 2:
            new_created_from_peaks.append(created_from_peaks[i])
            i += 1
        created_from_peaks = new_created_from_peaks
    return created_from_peaks


def get_tensor_boundaries(impure_peak, relevant_comp, iter_offset):
    """
    Returns a tuple containing the leftest and rightest peak/component boundaries.
    The boundaries are corrected by the iteration offset if the iterative
    PARAFAC algorithm is used.
    """
    created_from_peaks = get_comp_peaks(relevant_comp)
    lefts = []
    rights = []
    for peak in created_from_peaks:
        offset = get_offset_peak_to_comp(peak, relevant_comp)
        lefts.append(peak.left - offset)
        rights.append(peak.right - offset)

    left_boundaries = lefts + [impure_peak.left]
    right_boundaries = rights + [impure_peak.right]

    left = min(left_boundaries)
    right = max(right_boundaries)

    right = right + abs(iter_offset)  # always extend to the right

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


def get_comp_peak_data_list(relevant_comp, boundaries, iter_offset):
    """
    Returns a list of maximum-aligned peak data of the comp peaks of the
    relevant components. Moreover, returns the shape of the component list
    (how many slices per component in the list).
    """
    ze_peaks = []
    comp_tensor_shape = tuple()
    created_from_peaks = get_comp_peaks(relevant_comp)
    comp_tensor_shape = comp_tensor_shape + (len(created_from_peaks),)
    for peak in created_from_peaks:
        offset = get_offset_peak_to_comp(peak, relevant_comp)
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


def get_zero_ext_impure_peak_data(impure_peak, boundaries, iter_offset):
    """
    Returns zero-extended (to the boundaries) peak data of the impure peak.
    The location of the peak is corrected by the iteration offset in case
    the iterative PARAFAC algorithm is applied.
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


def normalize_peak_data(parafac_data):
    """
    Normalizes all slices to the maximum absorbance of the slice.
    """
    max_impure = sum_absorbance_by_time(parafac_data[-1]).max()
    data_normalized = []
    for data in parafac_data[:-1]:
        max_data = sum_absorbance_by_time(data).max()
        data_normalized.append(data / max_data * max_impure / 2)
    data_normalized.append(parafac_data[-1])
    return data_normalized


def create_data_tensor(parafac_data):
    """
    Takes a list of peak data and returns a data tensor in the format used for
    PARAFAC decomposition.
    """
    # add all data to data_flattened and create data tensor
    data_flattened = []
    for data in parafac_data:
        data_flattened.append(data[:, :, np.newaxis])
    data_tensor = np.concatenate(data_flattened, axis=2)
    return data_tensor


def get_parafac_tensor(impure_peak, quali_comp_db, iter_offset,
                       show_parafac_analytics):
    """
    Processor function of the data tensor creation. Takes in the impure peak,
    the qualitative component db and an iteration offset and returns a
    DataTensor object which is used for PARAFAC decomposition.
    """
    relevant_comp = get_relevant_comp(impure_peak, quali_comp_db)

    boundaries = get_tensor_boundaries(impure_peak, relevant_comp, iter_offset)

    comp_peaks, comp_tensor_shape = get_comp_peak_data_list(relevant_comp,
                                                            boundaries,
                                                            iter_offset)

    impure_peak, y_offset = get_zero_ext_impure_peak_data(impure_peak, boundaries,
                                                          iter_offset)

    parafac_data = comp_peaks + [impure_peak]
    normalized_data = normalize_peak_data(parafac_data)
    data_tensor = create_data_tensor(normalized_data)

    if show_parafac_analytics:
        print(f"new data tensor with boundaries {boundaries} and shape"
              f"{comp_tensor_shape} built from {relevant_comp}.")

    return DataTensor(data_tensor, boundaries, relevant_comp,
                      comp_tensor_shape, y_offset)

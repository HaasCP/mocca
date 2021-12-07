# flake8: noqa
from mocca.peak.models import BasePeak
from scipy.signal import find_peaks, peak_widths
import numpy as np


def trim_data(data, time, length):
    """Trims the 2D DADData in the time dimension to the length provided"""
    if length < data.shape[1]:
        return data[:, :length], time[:length]
    else:
        return data, time


def get_compound_names(compound_input):
    # iterate through all values to find compounds with nonzero concs
    # input is a dictionary mapping compounds to concentrations
    return [compound for compound, conc in compound_input.items() if conc != 0]


def absorbance_sum_by_time(data):
    """
    Sums the absorbances for each time point over all wavelengths

    Parameters
    ----------
    data : numpy.ndarray
        Actual experimental data with shape [# of wavelengths] x [timepoints].
        Generated from dataframe with absorbance_to_array function

    Returns
    -------
    numpy.ndarray
        A 1D array containing the sum of wavelengths at each time point
    """
    return data.sum(axis=0)


def get_peaks(data, absorbance_threshold, peaks_high_pass=None,
              peaks_low_pass=None, expand_peaks=True):
    """
    Finds all peaks of data.

    Parameters
    ----------
    data : numpy.ndarray
        Actual experimental data with shape [# of wavelengths] x [timepoints].
        Generated from dataframe with absorbance_to_array function

    absorbance_threshold : float
        The threshold below which peaks will. In other words, at at least one
        (wavelength, timepoint) will have absorbance greater than
        absorbance_threshold in order to be counted as a peak.

    peaks_high_pass : float
        TODO

    peaks_low_pass : float
        TODO

    expand_peaks : boolean
        If True, then peaks will be expanded to their peak boundaries. If this is
        set to False, then only timepoints with cumulative absorbance greater
        than absorbance_threshold will be counted as part of the peak.

    Returns
    -------
    peaks : list
        List of all peaks, as a list of tuples (maximum, left, right)
    """

    # TODO
    summed_data = absorbance_sum_by_time(data)
    new_data_thresh = summed_data.copy()
    new_data_thresh[new_data_thresh < absorbance_threshold] = 0  # summed absorbance threshold

    # get peaks
    peaks = get_peak_locs(new_data_thresh)
    peaks = merge_peaks(summed_data, peaks)
    peaks = expand_peaks(summed_data, peaks, abs_threshold=absorbance_threshold / 15)

    # TODO: finish function


def get_peak_locs(summed_data):
    """
    Finds all peaks of data.

    Parameters
    ----------
    summed_data : numpy.ndarray
        A 1D array representing the absorbances over time. Best used on
        data that already had data below threshold zeroed (see function
        filter_absorbance_by_threshold).

    Returns
    -------
    peaks : list
        List of all peaks, as a list of BasePeak classes
    """

    peak_maxima, _ = find_peaks(summed_data)
    peak_borders = peak_widths(summed_data, peak_maxima, rel_height=1)
    # use tolist to avoid problems during JSON export
    peaks = list(zip(peak_maxima.tolist(),
                     np.floor(peak_borders[2]).astype(int).tolist(),
                     np.ceil(peak_borders[3]).astype(int).tolist()))
    # peaks is now a list of tuples (maximum, left border, right border)
    return [BasePeak(maximum=maximum, left=left, right=right)
            for maximum, left, right in peaks]


def merge_peaks(summed_data, peaks):
    """
    Merges overlapping peaks in the data.

    Parameters
    ----------
    summed_data : numpy.ndarray
        A 1D array representing the absorbances over time. Best used on
        data that already had data below threshold zeroed (see function
        filter_absorbance_by_threshold).

    peaks : list
        List of all peaks as BasePeak objects

    Returns
    -------
    new_peaks : list
        List of all peaks in dictionary format with keys maximum, left, and right.
        Peaks that overlap are merged together into one BasePeak.
    """

    # TODO: check new peak functions
    new_peaks = set()
    for cur_peak in peaks:
        overlapping_peaks = []
        for other in peaks:
            if peak_overlap(cur_peak, other):  # TODO: fix function from peak class
                overlapping_peaks.append(other)

        merged_peak = [-1, np.inf, 0]
        for peak in overlapping_peaks:  # otherwise, merge peaks
            merged_peak[0] = peak['maximum'] if summed_data[peak['maximum']] > summed_data[merged_peak[0]] or merged_peak[0] == -1 else merged_peak[0]
            merged_peak[1] = min(merged_peak[1], peak['left'])
            merged_peak[2] = max(merged_peak[2], peak['right'])
        new_peaks.add(tuple(merged_peak))

    new_peaks = sorted(new_peaks, key=lambda x: x[0])
    peak_list = []
    for peak in new_peaks:
        peak_list.append(peak_to_dict(peak))  # TODO: fix

    if peak_list == peaks:
        return peak_list  # if merge does nothing, then return
    elif len(merge_peaks(summed_data, peak_list)) < len(peak_list):  # else, if the merging results in fewer peaks, recurse
        return merge_peaks(summed_data, peak_list)  # recursion takes care of some special cases of overlapping peaks, and merges the previous list if necessary
    return peak_list  # otherwise, just return fully merged list


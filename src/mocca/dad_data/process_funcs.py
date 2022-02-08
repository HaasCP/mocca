import numpy as np
from scipy.signal import find_peaks, peak_widths

from mocca.peak.models import BasePeak, PickedPeak
from mocca.chromatogram.model import Chromatogram
from mocca.chromatogram.utils import check_overlap

from mocca.dad_data.utils import sum_absorbance_by_time


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
    return [BasePeak(maximum=maximum, left=left, right=right, offset=0)
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

    new_peaks = set()
    for cur_peak in peaks:
        overlapping_peaks = []
        for other in peaks:
            if check_overlap(cur_peak, other):
                overlapping_peaks.append(other)

        merged_peak = [-1, np.inf, 0]
        for peak in overlapping_peaks:  # otherwise, merge peaks
            merged_peak[0] = (peak.maximum if summed_data[peak.maximum] >
                              summed_data[merged_peak[0]] or merged_peak[0] == -1
                              else merged_peak[0])
            merged_peak[1] = min(merged_peak[1], peak.left)
            merged_peak[2] = max(merged_peak[2], peak.right)
        new_peaks.add(tuple(merged_peak))

    new_peaks = sorted(new_peaks, key=lambda x: x[0])
    peak_list = [BasePeak(maximum=maximum, left=left, right=right, offset=0)
            for maximum, left, right in new_peaks]

    if peak_list == peaks:
        return peak_list  # if merge does nothing, then return
    elif len(merge_peaks(summed_data, peak_list)) < len(peak_list):  # else, if the merging results in fewer peaks, recurse
        return merge_peaks(summed_data, peak_list)  # recursion takes care of some special cases of overlapping peaks, and merges the previous list if necessary
    return peak_list  # otherwise, just return fully merged list


def pick_peaks(compound_data, experiment, absorbance_threshold,
               peaks_high_pass, peaks_low_pass):
    """
    Finds all peaks of data and returns them as a chromatogram

    Parameters
    ----------
    data : numpy.ndarray
        Actual experimental data with shape [# of wavelengths] x [timepoints].
        Generated from dataframe with absorbance_to_array function

    absorbance_threshold : float
        The threshold below which peaks will. In other words, at at least one
        (wavelength, timepoint) will have absorbance greater than
        absorbance_threshold in order to be counted as a peak.

    peaks_high_pass: float
        Time high pass filter only using peaks with a retention time
        greater than the here given value for data analysis
    
    peaks_low_pass : float
        Time low pass filter only using peaks with a retention time lower
        than the here given value for data analysis

    expand_peaks : boolean
        If True, then peaks will be expanded to their peak boundaries. If this is
        set to False, then only timepoints with cumulative absorbance greater
        than absorbance_threshold will be counted as part of the peak.

    Returns
    -------
    peaks : list
        List of all peaks, as a list of tuples (maximum, left, right)
    """

    summed_data = sum_absorbance_by_time(compound_data.data)
    new_data_thresh = summed_data.copy()
    new_data_thresh[new_data_thresh < absorbance_threshold] = 0  # summed absorbance threshold
    
    peak_locs = get_peak_locs(new_data_thresh)
    merged_peaks = merge_peaks(summed_data, peak_locs)
    
    if peaks_high_pass:
        merged_peaks = [peak for peak in merged_peaks if 
                        compound_data.time[peak.maximum] >= peaks_high_pass]
    if peaks_low_pass:
        merged_peaks = [peak for peak in merged_peaks if 
                        compound_data.time[peak.maximum] <= peaks_low_pass]
    
    merged_peaks = sorted(merged_peaks, key=lambda peak: peak.maximum)
    
    chromatogram = Chromatogram(experiment, compound_data)
    for idx, peak in enumerate(merged_peaks):
        chromatogram.insert_peak(PickedPeak(left= peak.left,
                                            right=peak.right,
                                            maximum=peak.maximum,
                                            offset=peak.offset,
                                            dataset=compound_data,
                                            idx=idx+1))
    return chromatogram

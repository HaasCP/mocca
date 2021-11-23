import numpy as np

from dataclasses import dataclass
import typing

from mocca.peak_purity import PeakPurityPredictor


@dataclass
class Peak():
    left : int
    right : int
    maximum : int
    dataset : typing.Any  # DADData parent of peak
    saturation : bool = None
    pure : bool = None
    istd : bool = None
    integral : float = None
    compound_id : str = None
    idx : int = None
    concentration : float = None

    def get_attribute(self, attribute):
        """
        Returns the attribute corresponding to the given string
        """
        return getattr(self, attribute)


class PeakProcessor():
    """
    Class which collects functions to process Peak objects.
    """
    def __init__(self, peak, absorbance_threshold, detector_limit,
                 print_purity_check=False):
        """
        Parameters
        ----------
        peak : mocca Peak object

        absorbance_threshold : int
            Minimum absorbance value when picking peaks

        detector_limit : int
            Absorbance values above which detector saturation is expected

        print_purity_check : bool, default=False
            Indicates if a print out of the peak purity predictor is given
        """
        if type(peak) == Peak:
            self.peak = peak
        else:
            raise Exception("The given peaks is not an object of the mocca Peak class")
        self.expand_threshold = absorbance_threshold / 20
        self.detector_limit = detector_limit
        self.print_purity_check = print_purity_check

    def expand_peak(self):
        """
        Expands peak boundaries to those actually in the data. It keeps expanding
        them until the absorbance falls below either abs_threshold or
        rel_threshold times the maximum overall absorbance sum.

        Modifies
        --------
        self.peak.left, self.peak.right : The left and right boundaries of the
        peak are updated
        """
        # sum absorbances over all wavelengths
        data = np.sum(self.peak.dataset.data, axis=0)
        # smoothing filter on time axis with a window length of 5
        data = np.convolve(data, np.ones(5), 'same') / 5

        prev_val = np.inf
        while data[self.peak.left] > self.expand_threshold and \
                prev_val > data[self.peak.left] and self.peak.left >= 0:
            prev_val = data[self.peak.left]
            self.peak.left -= 1

        if prev_val != np.inf:  # if peak was expanded, fix boundary, else don't change
            self.peak.left += 1

        prev_val = np.inf
        while data[self.peak.right] > self.peak.expand_threshold and \
                prev_val > data[self.peak.right] and self.peak.right <= len(data) - 1:
            prev_val = data[self.peak.right]
            self.peak.right += 1

        if prev_val != np.inf:  # if peak was expanded, fix boundary, else don't change
            self.peak.right -= 1

    def integrate_peak(self):
        """
        Integrates the peak.

        Modifies
        --------
        self.peak.integral : Integral of the peak
        """
        peak_data = self.peak.dataset.data[:, self.peak.left:self.peak.right]
        self.peak.integral = np.sum(peak_data).tolist()

    def check_peak_saturation(self):
        """
        Integrates the peak and sets self.integral to be that value.

        Modifies
        --------
        self.peak.saturation : Sets peak attribute to either True or False
            based on if the peak absorbance exceeds detector_limit.
        """
        max_absorbance = self.peak.dataset.data[:, self.peak.maximum].max()
        self.peak.saturation = max_absorbance > self.detector_limit

    def check_peak_purity(self):
        """
        Predicts purity of the peak.

        Modifies
        --------
        self.peak.pure : bool based on purity prediction.
        """
        purity_predictor = PeakPurityPredictor(self.peak)
        if self.print_purity_check:
            purity_predictor.show_analytics()
        self.peak.pure = purity_predictor.predict_peak_purity()

    def process_peak(self):
        self.expand_peak()
        self.integrate_peak()
        self.check_peak()
        return self.peak


class PeakPair():
    """
    Objects consisting of two peaks. Methods which are applied on two peaks.
    """
    def __init__(self, peak_1, peak_2):
        """
        Takes two Peak objects and checks for correct type
        """
        if type(peak_1) == Peak and type(peak_2) == Peak:
            self.peak_1 = peak_1
            self.peak_2 = peak_2
        else:
            raise Exception("At least one of the given peaks is not an \
                            object of the mocca Peak class")

    def check_same_dataset(self):
        """
        Throws Exception if the two peaks are not from the same dataset.
        """
        if self.peak_1.dataset != self.peak_2.dataset:
            raise Exception("Peaks are not from the same dataset, \
                            when comparing peak {} and {}!".format(self.peak_1.idx,
                            self.peak_2.idx))

    def check_overlap(self):
        """
        Returns True if self overlaps with the peak 'other', and False otherwise.

        Throws Exception if the two peaks are not from the same dataset.
        """
        self.check_same_dataset()
        return self.peak_1.left <= self.peak_2.left <= self.peak_1.right \
            or self.peak_2.left <= self.peak_1.left <= self.peak_2.right

    def distance_to(self, other):
        """
        Returns the distance from the maxima of self and the other peak.

        Throws Exception if the two peaks are not from the same dataset.
        """
        self.check_same_dataset(other)
        return abs(self.peak_1.maximum - self.peak_2.maximum)

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

    # Peak processing functions
    def expand_peak(self, absorbance_threshold):
        """
        Expands peak boundaries to those actually in the data. It keeps expanding
        them until the absorbance falls below either expand_threshold or
        rel_threshold times the maximum overall absorbance sum.

        Parameters
        ----------
        absorbance_threshold : int
            Minimum absorbance value when picking peaks

        Modifies
        --------
        self.peak.left, self.peak.right : The left and right boundaries of the
        peak are updated
        """
        expand_threshold = absorbance_threshold / 20
        # sum absorbances over all wavelengths
        data = np.sum(self.dataset.data, axis=0)
        # smoothing filter on time axis with a window length of 5
        data = np.convolve(data, np.ones(5), 'same') / 5

        prev_val = np.inf
        while data[self.left] > expand_threshold and \
                prev_val > data[self.left] and self.left >= 0:
            prev_val = data[self.left]
            self.left -= 1

        if prev_val != np.inf:  # if peak was expanded, fix boundary, else don't change
            self.left += 1

        prev_val = np.inf
        while data[self.right] > expand_threshold and \
                prev_val > data[self.right] and self.right <= len(data) - 1:
            prev_val = data[self.right]
            self.right += 1

        if prev_val != np.inf:  # if peak was expanded, fix boundary, else don't change
            self.right -= 1

    def integrate_peak(self):
        """
        Integrates the peak.

        Modifies
        --------
        self.peak.integral : Integral of the peak
        """
        peak_data = self.dataset.data[:, self.left:self.right]
        self.integral = np.sum(peak_data).tolist()

    def check_peak_saturation(self, detector_limit):
        """
        Integrates the peak and sets self.integral to be that value.
        Parameters
        ----------
        detector_limit : int
            Absorbance values above which detector saturation is expected

        Modifies
        --------
        self.peak.saturation : Sets peak attribute to either True or False
            based on if the peak absorbance exceeds detector_limit.
        """
        max_absorbance = self.dataset.data[:, self.maximum].max()
        self.saturation = max_absorbance > detector_limit

    def check_peak_purity(self, print_purity_check=False):
        """
        Predicts purity of the peak.
        Parameters
        ----------
        print_purity_check : bool, default=False
            Indicates if a print out of the peak purity predictor is given

        Modifies
        --------
        self.peak.pure : bool based on purity prediction.
        """
        purity_predictor = PeakPurityPredictor(self)
        if print_purity_check:
            purity_predictor.show_analytics()
        self.pure = purity_predictor.predict_peak_purity()

    def process_peak(self, absorbance_threshold, detector_limit,
                     print_purity_check=False):
        """
        Peak processing routine
        """
        self.expand_peak(absorbance_threshold)
        self.integrate_peak()
        self.check_peak_saturation(detector_limit)
        self.check_peak_purity(print_purity_check)

    # Functions including second peak
    def check_same_type(self, other):
        """
        Raises exception if given peak is not an mocca peak object
        """
        if type(other) != type(self):
            raise Exception("Given peak is not an \
                            object of the mocca Peak class")
        # TODO: handle exception and write test

    def check_same_dataset(self, other):
        """
        Raises Exception if the two peaks are not from the same dataset or same type.
        """
        self.check_same_type(other)
        if self.dataset != other.dataset:
            raise Exception("Peaks are not from the same dataset, \
                            when comparing peak {} and {}!".format(self.idx,
                            other.idx))
        # TODO: handle exception and write test

    def check_overlap(self, other):
        """
        Returns True if self overlaps with the peak 'other', and False otherwise.
        Raises Exception if the two peaks are not from the same dataset or same type.
        """
        self.check_same_dataset(other)
        return self.left <= other.left <= self.right \
            or other.left <= self.left <= other.right

    def distance_to(self, other):
        """
        Returns the distance from the maxima of self and the other peak.
        Raises Exception if the two peaks are not from the same dataset or same type.
        """
        self.check_same_dataset(other)
        return abs(self.maximum - other.maximum)

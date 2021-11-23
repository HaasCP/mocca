import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from mocca.utils import is_unimodal
from dataclasses import dataclass
import typing


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

# All following functions should be put in other class

    def check_same_dataset(self, other):
        """
        Throws Exception if the two peaks are not from the same dataset.
        """
        if self.dataset != other.dataset:
            raise Exception("Peaks are not from the same dataset, \
                            when comparing peak {} and {}!".format(self.idx, other.idx))

    def integrate_peak(self):
        """
        Integrates the peak and sets self.integral to be that value.
        """
        self.integral = np.sum(self.dataset.data[:, self.left:self.right]).tolist()

    def peak_overlap(self, other):
        """
        Returns True if self overlaps with the peak 'other', and False otherwise.

        Throws Exception if the two peaks are not from the same dataset.
        """
        self.check_same_dataset(other)
        return self.left <= other.left <= self.right \
            or other.left <= self.left <= other.right

    def distance_to(self, other):
        """
        Returns the distance from the maxima of self and the other peak.

        Throws Exception if the two peaks are not from the same dataset.
        """
        self.check_same_dataset(other)
        return abs(self.maximum - other.maximum)

    def check_peak(self, detector_limit, wavelength_filter=True,
                   data_filter=True, show_analytics=False):
        """
        Parameters
        ----------
        wavelength_filter : boolean, optional
            If true, filters out bad wavelengths (those with low absorbance and
            those exceeding the absorbance threshold). Default is true.

        data_filter : boolean, optional
            If true, filters out timepoints that are below the absorbance threshold.
            Default is true.

        show_analytics : boolean, optional
            If true, then prints plots and print statements for debugging. Default is
            false.

        Modifies
        --------
        self.pure
            Sets self.pure to either True or False based on prediction.

        self.saturation
            Sets self.saturation to either True or False based on if the peak
            absorbance exceeds detector_limit.

        """
        self.pure = self._check_peak_purity(detector_limit, wavelength_filter,
                                            data_filter, show_analytics)
        self.saturation = self._check_peak_saturation(detector_limit)

    def expand_peak(self, abs_threshold=30, rel_threshold=0.005):
        """
        Expands peak boundaries to those actually in the data. It keeps expanding
        them until the absorbance falls below either abs_threshold or
        rel_threshold times the maximum overall absorbance sum.

        Modifies
        --------
        self.left, self.right : The left and right boundaries of the peak are updated
        """
        data = np.sum(self.dataset.data, axis=0)
        data = np.convolve(data, np.ones(5), 'same') / 5
        # averaging filter of length 5

        # ensure peak boundaries are relatively consistent at different concentrations
        absorbance_threshold = min(abs_threshold, data[self.maximum] * rel_threshold)

        prev_val = np.inf
        while data[self.left] > absorbance_threshold and \
                prev_val > data[self.left] and self.left >= 0:
            prev_val = data[self.left]
            self.left -= 1

        if prev_val != np.inf:  # if peak was expanded, fix boundary, else don't change
            self.left += 1

        prev_val = np.inf
        while data[self.right] > absorbance_threshold and \
                prev_val > data[self.right] and self.right <= len(data) - 1:
            prev_val = data[self.right]
            self.right += 1

        if prev_val != np.inf:  # if peak was expanded, fix boundary, else don't change
            self.right -= 1

    def check_database(self, component_database, abs_threshold=100,
                       rel_threshold=0.01, similarity_threshold=0.9):
        """
        Computes if the peak matches one already seen in the database, based
        on retention time and UV-Vis spectra.

        Parameters
        ----------
        component_database : micdrop.ComponentDatabase
            A component database containing all seen components.

        abs_threshold : numeric
            The absolute number of timepoints that the retention time is allowed to
            be off by. Default is 100.

        rel_threshold : numeric
            The relative fraction of timepoints that the retention time is allowed to
            be off by. Default is 0.01.

        similarity_threshold : numeric
            The correlation needed between spectra to match a component in the
            component database. Default is 0.9

        Returns
        --------
        str : The name of the closest matching Component if there was a match,
        or "Unknown" if there was no match.
        """

        matches = []
        for component in component_database:
            similarity = np.corrcoef(component.spectra,
                                     self.dataset.data[:, self.maximum])[1, 0]
            distance_to_max = abs(component.maximum - self.maximum)
            if distance_to_max < abs_threshold + rel_threshold * component.maximum and \
                    similarity > similarity_threshold:
                matches.append({'name': component.name, 'similarity': similarity})
            # sort in decreasing order, based on similarity of spectra
            matches.sort(key=lambda x: -x['similarity'])

        if len(matches) == 0:
            print("Warning: check_database() found no matches for peak {}, setting compound_id as Unknown".format(self.idx))
            return "Unknown"

        if len(matches) > 1:
            print("Warning: check_database() found multiple matches for peak {}, the full list is {}".format(self.idx, matches))

        return matches[0]['name']

    def set_compound_id(self, component_database, abs_threshold=100,
                        rel_threshold=0.01, similarity_threshold=0.9):
        """
        Checks if the peak is in the database, and if so, sets self.compound_id.
        See check_database().

        Modifies
        --------
        self.compound_id
            If a match is found, then the peak's compound_id attribute is set.
        """

        if self.pure is None or not self.pure:
            print("Warning: Running set_compound_id() on impure peak {}.".format(self.idx))

        self.compound_id = self.check_database(component_database=component_database,
                                               abs_threshold=abs_threshold,
                                               rel_threshold=rel_threshold,
                                               similarity_threshold=similarity_threshold)

    def quantify_peak(self, quantification_database):
        """
        Computes the concentration of the compound in this peak.

        The attributes self.integral and self.compound_id must be set beforehand
        (through functions self.integrate_peak() and self.check_database())
        in order to quantify.

        Parameters
        ----------
        component_database : micdrop.component_database
            A component database containing all seen components.

        abs_threshold : numeric
            The absolute number of timepoints that the retention time is allowed to
            be off by. Default is 100.

        rel_threshold : numeric
            The relative fraction of timepoints that the retention time is allowed to
            be off by. Default is 0.01.

        Raises Exception if the attributes self.integral or self.compound_id are
        not set. Prints a text warning if self.pure is not set.
        """
        if self.integral is None or self.compound_id is None:
            raise Exception("Must set peak integral and compound_id before attempting to quantify!")
        raise NotImplementedError("TODO")

    def _check_peak_saturation(self, detector_limit):
        return self.dataset.data[:, self.maximum].max() > detector_limit

    def _filter_peak_data(self, data_filter, wavelength_filter, detector_limit):
        # helper for peak purity check
        peak_data = self.dataset.data[:, self.left:self.right]
        if wavelength_filter:
            wavelength_sums = peak_data.sum(axis=1)
            # filter out wavelengths whose max absorbance is < 0.01 max peak absorbance
            min_filter = wavelength_sums > np.max(wavelength_sums) * 0.01
            # filter out wavelengths whose max absorbance above saturation threshold
            max_filter = np.max(peak_data, axis=1) < detector_limit
            peak_data = peak_data[min_filter & max_filter]

        # only include timepoints whose absorbance sum is > 0.05 that of the maximum
        if data_filter:
            peak_data = peak_data[:, np.sum(peak_data, axis=0) > 0.05 * np.max(np.sum(peak_data, axis=0))]

        return peak_data

    def _check_peak_thresholds(self, peak_data, noise_variance, show_analytics):
        # helper for peak purity check
        peak_max_location = np.argmax(np.sum(peak_data, axis=0))
        correls = [(np.corrcoef(peak_data[:, i], peak_data[:, peak_max_location])[0, 1])**2
                   for i in range(peak_data.shape[1])]

        #  2.5 -> 0.5 in original agilent threshold
        agilent_thresholds = [(max(0, 1 - 2.5 * (noise_variance / np.var(peak_data[:, i]) +
                                                 noise_variance / np.var(peak_data[:, peak_max_location]))))**2
                                                 for i in range(peak_data.shape[1])]

        # check if > 90% of the points are greater than the modified agilent threshold.
        agilent_test = np.sum(np.greater(correls, agilent_thresholds)) / peak_data.shape[1]
        # if peak is pure, overall correlation across all relevant peaks should be high
        correls_test_1 = np.min(correls)
        # if peak is pure, average correlation across all peaks should be high
        correls_test_2 = np.mean(correls)
        # averaging filter of length 3
        # https://stackoverflow.com/questions/14313510/how-to-calculate-
        # rolling-moving-average-using-numpy-scipy
        unimodality_test = is_unimodal(np.convolve(correls, np.ones(3),
                                                   'valid') / 3, 0.99)

        pca = PCA(n_components=1)
        pca.fit(peak_data)
        #  if peak is pure, % explained with just one component should be high
        pca_test = pca.explained_variance_ratio_[0]

        if show_analytics:
            plt.plot(correls)
            plt.plot(agilent_thresholds)
            plt.show()
            for i in range(peak_data.shape[1]):
                plt.plot(peak_data[:, i])
            plt.show()
            print(f"Agilent Threshold (True for >0.9): {agilent_test} \n"
                  f"Unimodality Test (False for False): {unimodality_test} \n"
                  f"PCA Variance Explained (True for >0.995): {pca_test} \n"
                  f"Minimum Correlation (False for <0.9): {correls_test_1} \n"
                  f"Minimum Correlation (True for >0.95): {correls_test_1} \n"
                  f"Average Correlation (True for >0.98): {correls_test_2} \n")

        #  if agilent threshold reached, then probably pure
        if agilent_test > 0.9:
            return True
        #  for pure peak, correlation array emperically expected to be unimodal
        if not unimodality_test:
            return False
        #  if pca big enough, then probably pure
        if pca_test > 0.995:
            return True
        #  if any correlation is < 0.9, then probably impure somewhere
        if correls_test_1 < 0.9:
            return False
        #  otherwise, check if correlation shows that it is reasonably pure
        if correls_test_1 > 0.95:
            return True
        if correls_test_2 > 0.98:
            return True
        return False

    def _check_peak_purity(self, detector_limit, wavelength_filter,
                           data_filter, show_analytics):
        peak_data = self._filter_peak_data(data_filter=data_filter,
                                           wavelength_filter=wavelength_filter,
                                           detector_limit=detector_limit)
        #  filtered dataset with only timepoints whose max absorbance at
        #  any wavelength is 0.01 * max absorbance
        noise_data = self.dataset.data[:, np.max(self.dataset.data, axis=0) <
                                       0.01 * np.max(self.dataset.data)]
        #  take the average of the variance over all wavelengths
        noise_variance = np.mean(np.var(noise_data, axis=0))
        return self._check_peak_thresholds(peak_data=peak_data,
                                           noise_variance=noise_variance,
                                           show_analytics=show_analytics)

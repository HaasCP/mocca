import numpy as np
from sklearn.decomposition import PCA
from tensorly.decomposition import non_negative_parafac_hals
from mocca.peak.utils import get_peak_data
from mocca.peak.models import ParafacPeak, ProcessedPeak
from mocca.dad_data.models import ParafacData, DadData
from mocca.components.databases import QualiComponentDatabase
from mocca.components.models import QualiComponent
from mocca.components.utils import get_closest_peak
from typing import Optional, Tuple, List
import matplotlib.pyplot as plt

def get_zero_extended_peak_data(peak: ProcessedPeak, boundaries: Tuple[int, int]):
	"""
	Returns zero-extended peak data of shape (wavelengths, boundaries[1] - boundaries[0] + 1)
	"""
	peak_data = get_peak_data(peak)
	output_data = np.zeros((peak_data.shape[0], boundaries[1] - boundaries[0] + 1))
	output_data[:, peak.left - boundaries[0] : peak.right - boundaries[0] + 1] = peak_data
	return output_data


def get_relevant_data(left: int, right: int, component_db: QualiComponentDatabase) -> Tuple[List[np.ndarray], Tuple[int, int], int]:
	"""
	Given the left and right boundaries of the impure peak, scans the component database
	for any components that overlap with the impure peak, and chooses the peak most similar
	to average. Zero-extends the peak data until all peak data has the same length, and
	returns a list containing all relevant peaks, the overall peak boundaries, as well as
	the number of distinct peak components the data is from.
	"""

	def overlap(peak: ProcessedPeak, component: QualiComponent):
		return component.left <= peak.left <= component.right \
		    or peak.left <= component.left <= peak.right

	relevant_peaks = [peak for component in component_db if overlap(peak := get_closest_peak(component), component)]
	peak_boundaries = min(peak.left for peak in relevant_peaks), max(peak.right for peak in relevant_peaks)

	return [get_zero_extended_peak_data(peak=peak, boundaries=peak_boundaries) for peak in relevant_peaks], peak_boundaries, len(relevant_peaks)

def parafac(impure_peak: ProcessedPeak, component_db: QualiComponentDatabase) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Runs PARAFAC on an impure peak
    """

    # peak data is all the peaks that overlap with the impure peak in question
    peak_data, peak_boundaries, num_relevant_peaks = get_relevant_data(impure_peak.left, impure_peak.right, component_db)
    
    data_flattened = []

    # add all data to data_flattened and create data tensor
    for data in (get_zero_extended_peak_data(peak=impure_peak, boundaries=peak_boundaries) + peak_data):
        data_flattened.append(data[:, :, np.newaxis])

    data_tensor = np.concatenate(data_flattened, axis=2)
    ncomponents = max(num_relevant_peaks, estimate_num_peak_components(data_tensor.reshape(get_peak_data(impure_peak).shape[0], -1)))
    tensor = non_negative_parafac_hals(data_tensor, rank=ncomponents, init='svd', n_iter_max=1000, verbose = False, tol=1e-9)

    spectra = tensor[1][0]
    elutions = tensor[1][1]
    concs = tensor[1][2]

    spectral_normalization = np.sum(spectra, axis=0)
    elution_normalization = np.sum(elutions, axis=0)
    return spectra / spectral_normalization, elutions / elution_normalization, concs * elution_normalization * spectral_normalization


def estimate_num_components(data) -> int: 
    pca = PCA(n_components=10)
    peak_data_pca = pca.fit_transform(data)

    threshold = 0.995

    cum_sum = 0
    for idx in range(10):
        cum_sum += pca.explained_variance_ratio_[idx]
        if cum_sum > threshold:
            return idx + 1


def make_new_peak(old_peak: ProcessedPeak, spectra: np.ndarray, elution: np.ndarray, integral: float, component_num: int) -> ProcessedPeak:
	"""
	Makes a new ProcessedPeak corresponding to PARAFAC-processed impure processed peak
	Spectra, Elution, and Integral generated from PARAFAC

	If the old peak had index i, then the new peak has index i_{component_num}
	"""
	parafac_peak = ParafacPeak(left=old_peak.left,
						       right=old_peak.right,
						       maximum=old_peak.left + argmax(elution),
						       integral=integral,
						       spectra=spectra,
						       elution=elution,
						       offset=old_peak.offset)
	sim_parafac_data = ParafacData(parafac_peak=parafac_peak, original_dataset=old_peak.dataset)
	return ProcessedPeak(left=old_peak.left,
				         right=old_peak.right,
			    	     maximum=old_peak.left + argmax(elution),
			    	     dataset=sim_parafac_data,
			    	     idx=str(old_peak.idx) + '_' + component_num,
			    	     saturation=old_peak.saturation,
			    	     purity='PARAFAC',
					     integral=integral,
				         offset=old_peak.offset)

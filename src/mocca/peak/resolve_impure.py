import numpy as np
from sklearn.decomposition import PCA
from tensorly.decomposition import non_negative_parafac_hals
from mocca.peak.utils import get_peak_data
from mocca.peak.models import ParafacPeak, ProcessedPeak
from mocca.dad_data.models import ParafacData, DadData
from typing import Optional, Tuple
import matplotlib.pyplot as plt

def get_relevant_data(left: int, right: int, peak_database) -> list:
	raise NotImplementedError()
	pass


def estimate_num_components(data) -> int: 
    pca = PCA(n_components=10)
    peak_data_pca = pca.fit_transform(data)

    threshold = 0.995

    cum_sum = 0
    for idx in range(10):
        cum_sum += pca.explained_variance_ratio_[idx]
        if cum_sum > threshold:
            return idx + 1


def parafac_iter(input_peaks: list, left: int, right: int, ncomponents = Optional[int]) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Runs PARAFAC on a list of peaks
    Peak 0 is the actual impure peak in question
    Peak 1..n are the pure peaks
    """
    data_flattened = []
    for data in input_peaks:
        data_flattened.append(get_peak_data(peak)[:, :, np.newaxis])
        # TODO: how to run parafac with different peak left/right lengths?
        # TODO: what form should the data be in?

    data_tensor = np.concatenate(data_flattened, axis=2)
    ncomponents = max(2, estimate_num_peak_components(data_tensor.reshape(right - left + 1, -1)))
    tensor = non_negative_parafac_hals(data_tensor, rank=ncomponents, init='svd', n_iter_max=1000, verbose = False, tol=1e-9)

    spectra = tensor[1][0]
    elutions = tensor[1][1]
    concs = tensor[1][2]

    spectral_normalization = np.sum(spectra, axis=0)
    elution_normalization = np.sum(elutions, axis=0)
    return spectra / spectral_normalization, elutions / elution_normalization, concs * elution_normalization * spectral_normalization


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

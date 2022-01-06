import numpy as np

from sklearn.decomposition import PCA
from tensorly.decomposition import non_negative_parafac_hals
import matplotlib.pyplot as plt 

from mocca.peak.utils import get_peak_data
from mocca.peak.models import CorrectedPeak, IntegratedPeak

from mocca.dad_data.models import ParafacData


def check_comp_overlap(peak, comp):
    return comp.left <= peak.left <= comp.right \
        or peak.left <= comp.left <= peak.right


def get_relevant_comps(impure_peak, quali_comp_db):
    return [comp for comp in quali_comp_db if check_comp_overlap(impure_peak, comp)]


def get_compound_offset(created_from_peak, comp):
    return created_from_peak.maximum - comp.maximum


def get_parafac_boundaries(impure_peak, relevant_comps):
    lefts = []
    rights = []
    for comp in relevant_comps:
        created_from_peaks = comp.created_from
        for peak in created_from_peaks:
            offset = get_compound_offset(peak, comp)
            lefts.append(peak.left - offset)
            rights.append(peak.right - offset)

    left_boundaries = lefts + [impure_peak.left]
    right_boundaries = rights + [impure_peak.right]
    return (min(left_boundaries), max(right_boundaries))


def get_zeros_array(boundaries, n_wavelengths):
    return np.zeros((n_wavelengths, boundaries[1] - boundaries[0] + 1))


def get_zero_extended_peak_data(peak_data, left, boundaries):
    """
    Returns zero-extended peak data of shape (wavelengths, boundaries[1] - boundaries[0] + 1)
    """
    peak_data_ze = get_zeros_array(boundaries, peak_data.shape[0])
    rel_left = left - boundaries[0]
    rel_right = rel_left + peak_data.shape[1]
    peak_data_ze[:, rel_left:rel_right] = peak_data
    return peak_data_ze


def get_comp_ze_peaks(relevant_comps, boundaries):
    # max alignment!
    ze_peaks = []
    for comp in relevant_comps:
        created_from_peaks = comp.created_from
        for peak in created_from_peaks:
            offset = get_compound_offset(peak, comp)
            peak_data = get_peak_data(peak)
            peak_data_ze = get_zero_extended_peak_data(peak_data,
                                                       peak.left - offset,
                                                       boundaries)
            ze_peaks.append(peak_data_ze)
    return ze_peaks


def get_impure_ze_peak(impure_peak, boundaries):
    peak_data = get_peak_data(impure_peak)
    peak_data_ze = get_zero_extended_peak_data(peak_data, impure_peak.left,
                                               boundaries)
    return peak_data_ze

def get_parafac_data(impure_peak, quali_comp_db):
    
    relevant_comps = get_relevant_comps(impure_peak, quali_comp_db)
    
    boundaries = get_parafac_boundaries(impure_peak, relevant_comps)
    
    comp_peaks = get_comp_ze_peaks(relevant_comps, boundaries)
    
    impure_peak = get_impure_ze_peak(impure_peak, boundaries)
    
    return comp_peaks + [impure_peak], boundaries, relevant_comps


def estimate_num_components(data): 
    pca = PCA(n_components=10)
    _ = pca.fit_transform(data)

    threshold = 0.995

    cum_sum = 0
    for idx in range(10):
        cum_sum += pca.explained_variance_ratio_[idx]
        if cum_sum > threshold:
            return idx + 1


def parafac_analytics(normalized_spectra, normalized_elution, normalized_concs,
                      boundaries, relevant_comps):
    print(f"PARAFAC tensor has boundaries of {boundaries}.")
    plt.plot(normalized_spectra)
    for comp in relevant_comps:
        print(f"Compound used for PARAFAC data tensor: {comp.compound_id}")
        plt.plot([val / (max(comp.spectrum) / normalized_spectra.max()) 
                  for val in comp.spectrum], "--")
    plt.show()
    plt.plot(normalized_elution)
    plt.show()
    plt.plot(normalized_concs)
    plt.show()


def parafac(impure_peak, quali_comp_db, show_parafac_analytics=False):
    """
    Runs PARAFAC on an impure peak
    """

    parafac_data, boundaries, relevant_comps = get_parafac_data(impure_peak,
                                                                quali_comp_db)

    # add all data to data_flattened and create data tensor
    data_flattened = []
    for data in parafac_data:
        data_flattened.append(data[:, :, np.newaxis])
    data_tensor = np.concatenate(data_flattened, axis=2)
    
    pca_data = data_tensor.reshape(get_peak_data(impure_peak).shape[0], -1)
    pca_n_comps = estimate_num_components(pca_data)
    n_comps = max(len(relevant_comps), pca_n_comps)
    
    tensor = non_negative_parafac_hals(data_tensor, rank=n_comps,
                                       init='svd', n_iter_max=1000,
                                       verbose = False, tol=1e-9)

    spectra = tensor[1][0]
    elutions = tensor[1][1]
    integrals = tensor[1][2]

    spectral_normalization = np.sum(spectra, axis=0)
    elution_normalization = np.sum(elutions, axis=0)

    normalized_spectra = spectra / spectral_normalization
    normalized_elution = elutions / elution_normalization
    normalized_integrals = integrals * elution_normalization * spectral_normalization

    if show_parafac_analytics:
        parafac_analytics(normalized_spectra, normalized_elution,
                          normalized_integrals, boundaries, relevant_comps)

    parafac_tensor = (normalized_spectra, normalized_elution, normalized_integrals)
    
    return parafac_tensor, boundaries


def create_parafac_peaks(impure_peak, parafac_tensor, boundaries,
                         absorbance_threshold):
    """
    Makes a new ProcessedPeak corresponding to PARAFAC-processed impure processed peak
    Spectra, Elution, and Integral generated from PARAFAC

    If the old peak had index i, then the new peak has index i_{component_num}
    """
    n_comps = parafac_tensor[0].shape[1]
    if any(dim.shape[1] != n_comps for dim in parafac_tensor):
        raise ValueError("All dimensions of PARAFAC tensor must contain same "
                         "amount of components")

    parafac_peaks = []
    for i in range(n_comps):
        #  get tensor for one parafac comonent
        parafac_comp_tensor = (parafac_tensor[0][:, i],
                               parafac_tensor[1][:, i],
                               parafac_tensor[2][:, i])
        if type(impure_peak) == CorrectedPeak:
            parafac_peak = CorrectedPeak(left=boundaries[0],
                                         right=boundaries[1],
                                         maximum=(boundaries[0] +
                                                  np.argmax(parafac_comp_tensor[1])),
                                         dataset=ParafacData(impure_peak,
                                                             parafac_comp_tensor,
                                                             boundaries),
                                         idx=-impure_peak.idx,
                                         saturation=impure_peak.saturation,
                                         pure=True,
                                         integral=parafac_comp_tensor[2][-1], # reaction run is last in run dimension
                                         offset=0,
                                         istd=impure_peak.istd)
        elif type(impure_peak) == IntegratedPeak:
            parafac_peak = IntegratedPeak(left=boundaries[0],
                                         right=boundaries[1],
                                         maximum=(boundaries[0] +
                                                  np.argmax(parafac_comp_tensor[1])),
                                         dataset=ParafacData(impure_peak,
                                                             parafac_comp_tensor,
                                                             boundaries),
                                         idx=-impure_peak.idx,
                                         saturation=impure_peak.saturation,
                                         pure=True,
                                         integral=parafac_comp_tensor[2][-1])
        else:
            raise TypeError(f"Given impure peak is of type {type(impure_peak)}. "
                            "Only mocca IntegratedPeak and CorrectedPeak types "
                            "are allowed.")
        if np.max(parafac_peak.dataset.data) > absorbance_threshold:
            parafac_peaks.append(parafac_peak)
    return parafac_peaks


def get_parafac_peaks(impure_peak, quali_comp_db, absorbance_threshold,
                      show_parafac_analytics):
    parafac_tensor, boundaries = parafac(impure_peak, quali_comp_db,
                                         show_parafac_analytics=False)
    parafac_peaks = create_parafac_peaks(impure_peak, parafac_tensor, boundaries,
                                         absorbance_threshold)
    return parafac_peaks

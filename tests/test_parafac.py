from chromatogram_gen import generate_test_chromatograms, plot_test_data, test_spectra, generate_elution_profile
from mocca.peak.models import ParafacPeak, PickedPeak
from mocca.dad_data.process_funcs import pick_peaks
from mocca.dad_data.models import ParafacData
from mocca.peak.check import check_peak_purity
import matplotlib.pyplot as plt
import numpy as np

test_data = generate_test_chromatograms()
# plot_test_data(test_data[0])

chromatogram = pick_peaks(compound_data=test_data[0], absorbance_threshold=0.01, peaks_high_pass=None, peaks_low_pass=None)

def test_parafac_data_creation_1():
	peak = PickedPeak(left=750, right=899, maximum=850, dataset=test_data[0], idx=1)
	check_peak_purity(peak, show_analytics=False)
	spec1 = test_spectra[8] #get overlapping peak 1 spectra
	spec2 = test_spectra[9] #get overlapping peak 2 spectra
	elution1 = generate_elution_profile(elution_time=50, elution_width=15, length=150)
	elution2 = generate_elution_profile(elution_time=100, elution_width=15, length=150)

	# integral, spectra, elution normally calculated from parafac procedure
	peak1 = ParafacPeak(left=750, right=899, maximum=800, integral=1, spectra=spec1, elution=elution1, offset=0)
	peak2 = ParafacPeak(left=750, right=899, maximum=850, integral=1, spectra=spec2, elution=elution2, offset=0)

	parafac_data1 = ParafacData(hplc_system_tag=None, path=None, parafac_peak=peak1, original_dataset=test_data[0])
	parafac_data2 = ParafacData(hplc_system_tag=None, path=None, parafac_peak=peak2, original_dataset=test_data[0])

	# show peaks:
	# plt.plot(np.sum(test_data[0].data, axis=0))
	# plt.plot(np.sum(parafac_data1.data, axis=0))
	# plt.plot(np.sum(parafac_data2.data, axis=0))
	# plt.show()

	# sum of two deconvolved peaks should be approximately the peak in the original dataset
	assert abs(np.sum((parafac_data1.data + parafac_data2.data - test_data[0].data)[750:900])) < 0.05

# tests for peak.py
from mocca.peak import Peak
from mocca.component_database import ComponentDatabase
from mocca.quant_database import QuantificationDatabase
import matplotlib.pyplot as plt
import numpy as np
import pytest
from sim_test_data import test_spectra
from chromatogram_gen import Chromatogram
### TODO: make dataset actually point to an object of DADData class, rather than Chromatogram class

#### MODIFY AS NEEDED FOR TESTING
show_peak_purity_analytics = False


#### TEST DATA GENERATION
test_data_1 = Chromatogram(1000, concs=[1] * 10, elution_widths=[5, 10, 5, 5, 3, 10, 3, 6, 15, 15], elution_times = [100, 150, 200, 250, 300, 400, 500, 700, 800, 850])
test_data_2 = Chromatogram(1000) #random chromatogram
test_data_3 = Chromatogram(1000, concs=[100000] + [0] * 9, elution_times=[200] + [0]*9) #random chromatogram
test_data_4 = Chromatogram(1000, concs=[1, 0, 0, 0, 1, 0, 0, 0, 0, 0], elution_widths=[10]*10, elution_times=[50, 0, 0, 0, 55, 0, 0, 0, 0, 0]) #overlapping
test_data_5 = Chromatogram(1000, concs=[1, 0, 0, 0, 0.2, 0, 0.2, 0, 1, 0], elution_widths=[15]*10, elution_times=[50, 0, 0, 0, 51, 0, 53, 0, 60, 0]) #overlapping
test_data_6 = Chromatogram(1000, concs=[10000, 0, 0, 0, 2000, 0, 2000, 0, 10000, 0], elution_widths=[15]*10, elution_times=[50, 0, 0, 0, 51, 0, 53, 0, 60, 0]) #overlapping, saturated
test_data_7 = Chromatogram(1000, concs=[1] * 10, elution_widths=[5, 10, 5, 5, 3, 10, 3, 6, 15, 15], elution_times = [100, 150, 200, 250, 300, 400, 500, 700, 800, 850], add_noise=False)
test_data_8 = Chromatogram(1000, concs=[1,0,0,0,0,0,0,0,0,0], elution_widths=[10]*10, elution_times = [150]*10)


### TEST DATA VISUALIZATION
def plot_test_data(data):
    """
    data is an input of class Chromatogram

    For use in showing test data plots only, not for actual testing.
    """
    plt.plot(np.sum(data.data, axis = 0))
    plt.show()
    plt.plot(data.data)
    plt.show()

"""
plot_test_data(test_data_1)
plot_test_data(test_data_2)
plot_test_data(test_data_3)
plot_test_data(test_data_4)
plot_test_data(test_data_5)
"""

### ACTUAL TESTS

### PEAK CLASS TESTS
def test_peak_overlap_1():
    #check that non-overlapping peaks are detected as non-overlapping
    peak_1 = Peak(left = 100, right = 200, maximum = 150, dataset = test_data_1)
    peak_2 = Peak(left = 250, right = 350, maximum = 300, dataset = test_data_1)
    assert peak_1.peak_overlap(peak_2) == False
    assert peak_2.peak_overlap(peak_1) == False

def test_peak_overlap_2():
    #check that overlapping peaks are detected as overlapping
    peak_1 = Peak(left = 100, right = 200, maximum = 150, dataset = test_data_1)
    peak_2 = Peak(left = 150, right = 250, maximum = 200, dataset = test_data_1)
    assert peak_1.peak_overlap(peak_2) == True
    assert peak_2.peak_overlap(peak_1) == True

def test_peak_overlap_3():
    #check that different dataset peaks do throw an exception
    with pytest.raises(Exception):
        peak_1 = Peak(left = 100, right = 200, maximum = 150, dataset = test_data_1)
        peak_2 = Peak(left = 150, right = 250, maximum = 200, dataset = test_data_2)
        peak_1.peak_overlap(peak_2)

def test_peak_distance_1():
    #check that peak distances are correct
    peak_1 = Peak(left = 100, right = 200, maximum = 150, dataset = test_data_1)
    peak_2 = Peak(left = 250, right = 350, maximum = 300, dataset = test_data_1)
    assert peak_1.distance_to(peak_2) == 150
    assert peak_2.distance_to(peak_1) == 150

def test_peak_distance_2():
    #check that peak distances are correct
    with pytest.raises(Exception):
        peak_1 = Peak(left = 100, right = 200, maximum = 150, dataset = test_data_1)
        peak_2 = Peak(left = 250, right = 350, maximum = 300, dataset = test_data_2)
        peak_1.distance_to(peak_2)

def test_peak_integral_1():
    #check that peak integral is approximately 0 over empty peak area
    peak = Peak(left = 550, right = 650, maximum = 600, dataset = test_data_1)
    peak.integrate_peak()
    assert peak.integral is not None
    assert peak.integral < 0.005 #approximate noise level

def test_peak_integral_2():
    #check that peak integral is large over actual peak
    peak = Peak(left = 90, right = 100, maximum = 110, dataset = test_data_1)
    peak.integrate_peak()
    assert peak.integral is not None
    assert peak.integral > 0.9 #each peak should integrate to roughly its concentration in the chromatogram_gen (as spectra are normalized to area 1)

def test_peak_purity_saturation_1():
    peak = Peak(left = 50, right = 100, maximum = 150, dataset = test_data_1)
    peak.check_peak(detector_limit = 10, show_analytics = show_peak_purity_analytics, wavelength_filter = False, data_filter = False)
    assert peak.pure == True
    assert peak.saturation == False

def test_peak_purity_saturation_2():
    peak = Peak(left = 50, right = 100, maximum = 150, dataset = test_data_1)
    peak.check_peak(detector_limit = 10, show_analytics = show_peak_purity_analytics, wavelength_filter = True, data_filter = False)
    assert peak.pure == True
    assert peak.saturation == False

def test_peak_purity_saturation_3():
    peak = Peak(left = 50, right = 100, maximum = 150, dataset = test_data_1)
    peak.check_peak(detector_limit = 10, show_analytics = show_peak_purity_analytics, wavelength_filter = False, data_filter = True)
    assert peak.pure == True
    assert peak.saturation == False

def test_peak_purity_saturation_4():
    peak = Peak(left = 50, right = 100, maximum = 150, dataset = test_data_1)
    peak.check_peak(detector_limit = 10, show_analytics = show_peak_purity_analytics, wavelength_filter = True, data_filter = True)
    assert peak.pure == True
    assert peak.saturation == False

def test_peak_purity_saturation_5():
	#check for saturated peak
    peak = Peak(left = 150, right = 250, maximum = 200, dataset = test_data_3)
    peak.check_peak(detector_limit = 10, show_analytics = show_peak_purity_analytics)
    assert peak.pure == True
    assert peak.saturation == True

def test_peak_purity_saturation_6():
	#check for impure peak
    peak = Peak(left = 40, right = 65, maximum = 53, dataset = test_data_4)
    peak.check_peak(detector_limit = 10, show_analytics = show_peak_purity_analytics)
    assert peak.pure == False
    assert peak.saturation == False

def test_peak_purity_saturation_7():
	#check for impure peak
    peak = Peak(left = 40, right = 70, maximum = 55, dataset = test_data_5)
    peak.check_peak(detector_limit = 10, show_analytics = show_peak_purity_analytics)
    assert peak.pure == False
    assert peak.saturation == False

def test_peak_purity_saturation_8():
	#check for impure and saturated peak
    peak = Peak(left = 40, right = 70, maximum = 55, dataset = test_data_6)
    peak.check_peak(detector_limit = 10, show_analytics = show_peak_purity_analytics)
    assert peak.pure == False
    assert peak.saturation == True

def test_peak_expand_1():
	#expand big peak
    peak = Peak(left = 848, right = 852, maximum = 850, dataset = test_data_1)
    peak.expand_peak()
    assert peak.left < 835
    assert peak.right > 865

def test_peak_expand_2():
	#expand small peak
    peak = Peak(left = 298, right = 302, maximum = 300, dataset = test_data_1)
    peak.expand_peak()
    assert peak.left < 295
    assert peak.right > 305

def test_peak_expand_3():
	#does not expand correct peak
    peak = Peak(left = 285, right = 315, maximum = 300, dataset = test_data_1)
    peak.expand_peak()
    assert abs(peak.left - 285) <= 1 # +/- 1 error for tolerance due to noise
    assert abs(peak.right - 315) <= 1

def test_peak_expand_4():
	#does not expand correct peak, no noise dataset
    peak = Peak(left = 285, right = 315, maximum = 300, dataset = test_data_7)
    peak.expand_peak()
    assert peak.left == 285 
    assert peak.right == 315

### COMPONENT DATABASE CLASS TESTS

def test_component_database_1():
    #try adding peaks to component database
    peak_1 = Peak(left = 130, right = 170, maximum = 150, dataset = test_data_1)
    peak_2 = Peak(left = 280, right = 330, maximum = 300, dataset = test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')
    component_database.add_peak(peak_2, 'Component 2')
    assert 'Component 1' in component_database
    assert 'Component 2' in component_database
    assert 'Component 3' not in component_database
    assert 'Component' not in component_database
    assert component_database['Component 1'].left == 130
    assert component_database['Component 1'].right == 170
    assert component_database['Component 1'].maximum == 150
    assert np.corrcoef(component_database['Component 1'].spectra, test_spectra[1])[1,0] > 0.999 #check that spectra are same
    components = []
    for component in component_database:
        components.append(component)
    assert len(components) == 2

def test_component_database_2():
    #component database add peak with same name
    with pytest.raises(Exception):
        peak_1 = Peak(left = 130, right = 170, maximum = 150, dataset = test_data_1)
        peak_2 = Peak(left = 280, right = 330, maximum = 300, dataset = test_data_1)
        component_database = ComponentDatabase()
        component_database.add_peak(peak_1, 'Component 1')
        component_database.add_peak(peak_2, 'Component 1')

def test_component_database_3():
    #component database get unadded item
    with pytest.raises(Exception):
        peak_1 = Peak(left = 130, right = 170, maximum = 150, dataset = test_data_1)
        peak_2 = Peak(left = 280, right = 330, maximum = 300, dataset = test_data_1)
        component_database = ComponentDatabase()
        component_database.add_peak(peak_1, 'Component 1')
        component_database.add_peak(peak_2, 'Component 2')
        component_database['Component 3']

### PEAK CHECK DATABASE TESTS 

def test_check_database_1():
    #normal functioning database
    peak_1 = Peak(left = 130, right = 170, maximum = 150, dataset = test_data_1)
    peak_2 = Peak(left = 280, right = 330, maximum = 300, dataset = test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')
    component_database.add_peak(peak_2, 'Component 2')

    test_peak_1 = Peak(left = 120, right = 165, maximum = 152, dataset = test_data_7) #minor variations in peak location
    test_peak_2 = Peak(left = 270, right = 325, maximum = 298, dataset = test_data_7)

    test_peaks = [test_peak_1, test_peak_2]
    for peak in test_peaks:
        peak.check_peak(detector_limit = 10)
        peak.set_compound_id(component_database)
    
    assert test_peak_1.compound_id == 'Component 1'
    assert test_peak_2.compound_id == 'Component 2'

def test_check_database_2():
    #add different peaks at same retention time
    peak_1 = Peak(left = 130, right = 170, maximum = 150, dataset = test_data_8) #similar peak at same time point
    peak_2 = Peak(left = 130, right = 170, maximum = 150, dataset = test_data_1) 
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')
    component_database.add_peak(peak_2, 'Component 2')

    test_peak_1 = Peak(left = 120, right = 165, maximum = 152, dataset = test_data_7) #should match to peak 2
    test_peak_2 = Peak(left = 270, right = 325, maximum = 298, dataset = test_data_7) #should match to nothing

    test_peaks = [test_peak_1, test_peak_2]
    for peak in test_peaks:
        peak.check_peak(detector_limit = 10)
        peak.set_compound_id(component_database, similarity_threshold = 0.85) #custom similarity threshold to actually cause two matches
    
    assert test_peak_1.check_database(component_database) == 'Component 2'
    assert test_peak_2.check_database(component_database) == 'Unknown'
    assert test_peak_1.compound_id == 'Component 2'
    assert test_peak_2.compound_id == 'Unknown'

def test_peak_quantification_1():
    peak_1 = Peak(left = 130, right = 170, maximum = 150, dataset = test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')

    quantification_database = QuantificationDatabase()

    #error because self.integral and self.compound_id is None
    with pytest.raises(Exception):
        peak_1.quantify_peak(quantification_database)

    peak_1.integrate_peak()
    peak_1.set_compound_id(component_database)

    peak_1.quantify_peak(quantification_database)
    

# flake8: noqa


# All following functions should be put in other class





def test_peak_integral_1():
    # check that peak integral is approximately 0 over empty peak area
    peak = Peak(left=550, right=650, maximum=600, dataset=test_data[0])
    peak.integrate_peak()
    assert peak.integral is not None
    assert peak.integral < 0.005  # approximate noise level


def test_peak_integral_2():
    # check that peak integral is large over actual peak
    peak = Peak(left=90, right=110, maximum=100, dataset=test_data[0])
    peak.integrate_peak()
    assert peak.integral is not None
    assert peak.integral > 0.9
    # each peak should integrate to roughly its concentration in the chromatogram_gen
    # as spectra are normalized to area 1

# PEAK CHECK DATABASE TESTS

def test_check_database_1():
    # normal functioning database
    peak_1 = Peak(left=130, right=170, maximum=150, dataset=test_data_1)
    peak_2 = Peak(left=280, right=330, maximum=300, dataset=test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')
    component_database.add_peak(peak_2, 'Component 2')

    # minor variations in peak location
    test_peak_1 = Peak(left=120, right=165, maximum=152, dataset=test_data_7)
    test_peak_2 = Peak(left=270, right=325, maximum=298, dataset=test_data_7)

    test_peaks = [test_peak_1, test_peak_2]
    for peak in test_peaks:
        peak.process_peak(absorbance_threshold=0.1, detector_limit=10)
        peak.set_compound_id(component_database)

    assert test_peak_1.compound_id == 'Component 1'
    assert test_peak_2.compound_id == 'Component 2'


def test_check_database_2():
    # add different peaks at same retention time
    # similar peak at same time point
    peak_1 = Peak(left=130, right=170, maximum=150, dataset=test_data_8)
    peak_2 = Peak(left=130, right=170, maximum=150, dataset=test_data_1)
    component_database = ComponentDatabase()
    component_database.add_peak(peak_1, 'Component 1')
    component_database.add_peak(peak_2, 'Component 2')

    test_peak_1 = Peak(left=120, right=165,
                       maximum=152, dataset=test_data_7)  # peak 2
    test_peak_2 = Peak(left=270, right=325,
                       maximum=298, dataset=test_data_7)  # nothing

    test_peaks = [test_peak_1, test_peak_2]
    for peak in test_peaks:
        peak.process_peak(absorbance_threshold=0.1, detector_limit=10)
        peak.set_compound_id(component_database, similarity_threshold=0.85)
        # custom similarity threshold to actually cause two matches

    assert test_peak_1.check_database(component_database) == 'Component 2'
    assert test_peak_2.check_database(component_database) == 'Unknown'
    assert test_peak_1.compound_id == 'Component 2'
    assert test_peak_2.compound_id == 'Unknown'
# flake8: noqa


# All following functions should be put in other class

def get_compound_id(peak, component_db, spectrum_correl_coef_thresh,
                    relative_distance_thresh, print_out):
    """
    Checks if the peak is in the database, and if so, sets self.compound_id.
    See check_database().

    Modifies
    --------
    self.compound_id
        If a match is found, then the peak's compound_id attribute is set.
    """
    if peak.pure is True:
        predicted_compound_id = get_filtered_similarity_dicts(peak, component_db, 
                                                    spectrum_correl_coef_thresh,
                                                    relative_distance_thresh, 
                                                    print_out)
        if predicted_compound_id is None:
            component_db.increment_unkown_counter()
            return "unknown_" + str(component_db.unknown_counter)
        else:
            return predicted_compound_id
    else:
        return None

class QuantiComponentDatabase():
    def __init__(self):
        self.database = {}

    def add_compound_concentration(self, compound_id, conc_info):
        """
        Adds tuple conc_info = (integral, concentration) to the component
        database under compound name name. Creates a new dictionary entry
        if that compound was originally not present.
        """
        if compound_id not in self:
            self.database[compound_id] = []
        self.database[compound_id].append((conc_info[0], conc_info[1]))

    def __contains__(self, compound_id):
        """Allows statements such as "'A' in quantification_database"
        to see if that Component is inside the quantification database"""
        return compound_id in self.database

    def quantify_peak(self, integral, compound_id):
        """
        Quantifies specified peak, given peak integral and compound id.
        """
        if compound_id not in self:
            raise Exception("Compound {} not found in Quantification Database!".format(compound_id))
        integrals, concentrations = zip(*self.database[compound_id])
        return np.linalg.lstsq(np.reshape(integrals, (-1, 1)), concentrations, rcond=0)[0][0] * integral
    
    def quantify_peak(self, quantification_database):
        """
        Computes the concentration of the compound in this peak.

        The attributes self.integral and self.compound_id must be set beforehand
        (through functions self.integrate_peak() and self.check_database())
        in order to quantify.

        Parameters
        ----------
        quantification_database : micdrop.QuantificationDatabase
            A quantification database containing all seen components.

        Raises Exception if the attributes self.integral or self.compound_id are
        not set. Prints a text warning if self.pure is not set.

        Modifies
        --------
        self.concentration : sets concentration to that predicted by integral
        """
        if self.integral is None or self.compound_id is None:
            raise Exception("Must set peak integral and compound_id before \
                            attempting to quantify!")

        if self.pure is None or not self.pure:
            logging.warning("Warning: Running quantify_peak() on impure peak \
                            {}.".format(self.idx))

        self.concentration = quantification_database.quantify_peak(self.integral,
                                                                   self.compound_id)


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
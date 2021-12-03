# flake8: noqa
import numpy as np


class QuantificationDatabase():
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

from dataclasses import dataclass


@dataclass
class Component():
    compound_id: str
    left: int
    right: int
    maximum: int
    spectra: list


class ComponentDatabase():
    def __init__(self):
        self.database = {}

    def add_peak(self, peak, compound_id):
        """Creates a Component with name compound_id and data from the
        peak (of Peak class) and adds it into the database"""
        if compound_id in self:
            raise AttributeError("Component Database already contains compound {}".format(compound_id))
        peak = Component(compound_id=compound_id, left=peak.left, right=peak.right,
                         maximum=peak.maximum,
                         spectra=peak.dataset.data[:, peak.maximum].tolist())
        self.database[compound_id] = peak

    def __contains__(self, compound_id):
        """Allows statements such as "'Component 1' in component_database"
        to see if that Component is inside the database"""
        return compound_id in self.database

    def __getitem__(self, compound_id):
        """Allows statements such as "component_database['Component 1']" to access
           that Component"""
        if compound_id in self:
            return self.database[compound_id]
        else:
            raise AttributeError("{} not found in Component Database!".format(compound_id))

    def __iter__(self):
        """Yields all components inside the database.
           Allows statements such as "for component in component_database:" """
        for component in self.database.values():
            yield component

from dataclasses import dataclass

@dataclass
class Component():
    name: str
    left: int
    right: int
    maximum: int
    spectra: list

class ComponentDatabase():
    def __init__(self):
        self.database = {}

    def add_peak(self, peak, name):
        """Creates a Component with name name and data from the peak (of Peak class), and adds it into the database"""
        if name in self:
            raise AttributeError("Component Database already contains compound {}".format(name))
        peak = Component(name = name, left = peak.left, right = peak.right, maximum = peak.maximum, spectra = peak.dataset.data[:, peak.maximum].tolist())
        self.database[name] = peak

    def __contains__(self, name):
        """Allows statements such as 'Peak 1' in component_database to see if that Component is inside the database"""
        return name in self.database

    def __getitem__(self, name):
        """Allows statements such as component_database['Peak 1'] to access that Component"""
        if name in self:
            return self.database[name]
        else:
            raise AttributeError("{} not found in Component Database!".format(name))

    def __iter__(self):
        """Yields all components inside the database. Allows statements such as for component in component_database:"""
        for component in self.database.values():
            yield component
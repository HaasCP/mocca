from dataclasses import dataclass
from typing import List

from mocca.peak.models import ProcessedPeak


@dataclass
class QualiComponent():
    """
    Class of a qualitative component created from a number of peaks of the same
    compound_id.
    """
    compound_id : str
    left : int
    right : int
    maximum : int
    offset : int
    spectrum : list
    created_from : List[ProcessedPeak]

    def __repr__(self):
        kws = [f"{key}={value!r}" if key != "spectrum" else
               f"{key}={type(value)!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))


@dataclass
class QuantComponent():
    """
    Class of a quantitative component created from a number of peaks of the same
    compound_id.
    """
    compound_id : str
    calib_factors : dict
    calib_data : dict
    calib_scores : dict
    created_from : List[ProcessedPeak]

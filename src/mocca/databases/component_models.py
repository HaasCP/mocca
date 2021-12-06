from dataclasses import dataclass
from typing import List

from mocca.peak.models import ProcessedPeak


@dataclass
class QualiComponent():
    compound_id: str
    left: int
    right: int
    maximum: int
    spectrum: list
    created_from: List[ProcessedPeak]

"""
@dataclass
class QuantComponent():
    compound_id: str
    spectrum: list
    calibration_factor: float
    compound_peaks: List[ProcessedPeak]
    istd_peaks: List[ProcessedPeak] = []
"""
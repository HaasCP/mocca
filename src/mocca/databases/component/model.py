from dataclasses import dataclass
from typing import List

from mocca.peak.models import ProcessedPeak


@dataclass
class Component():
    compound_id: str
    left: int
    right: int
    maximum: int
    spectrum: list
    created_from: List[ProcessedPeak]

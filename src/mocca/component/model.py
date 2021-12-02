from dataclasses import dataclass
import logging
import numpy as np


@dataclass
class Component():
    compound_id: str
    left: int
    right: int
    maximum: int
    spectrum: list
    created_from: list

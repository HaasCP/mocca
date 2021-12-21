#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 09:27:43 2021

@author: haascp
"""

import os
from dataclasses import dataclass
from typing import Optional, List

from mocca.campaign.compound import Compound
from mocca.campaign.istd import InternalStandard

@dataclass()
class Experiment():
    """
    Data container to store user input.
    """
    path : str
    compound: Optional[Compound]
    istd: Optional[List[InternalStandard]]
    
    def __post_init__(self):
        if self.istd is not None and type(self.istd) != list:
            self.istd = [self.istd]
        if not os.path.exists(self.path):
            raise ValueError("Given path {} does not exist.".format(self.path))
        if self.compound.solvent and self.istd:
            raise ValueError("Solvent run has an internal standard added. Use "
                             "solvent == True only for pure solvent runs. These "
                             "runs will be analyzed first and should cover the "
                             "case that all samples are recorded in a UV-Vis "
                             "active solvent. Solvents can also be added as "
                             "compounds later on with solvent == False.")
        #solvents can be added also later via usual compound adding
        if self.compound.istd and self.istd:
            if self.compound.key == self.istd.key:
                raise ValueError("Internal standard cannot be analyzed relative "
                                 "to itself. If the internal standard should be "
                                 "added as compound, do not give internal "
                                 "standard parameter. If a run containing "
                                 "internal standard should be analyzed do not "
                                 "give internal standard as a compound.")

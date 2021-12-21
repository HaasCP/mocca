#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 12:00:20 2021

@author: haascp
"""

from dataclasses import dataclass
from typing import Optional

@dataclass()
class Compound():
    """
    Data container to store user input.
    """
    key: str
    conc: Optional[float] = None
    # following are only for ordering batch of runs
    solvent: bool = False
    istd: bool = False

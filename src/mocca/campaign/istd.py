#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 12:05:38 2021

@author: haascp
"""

from dataclasses import dataclass
from typing import Optional

@dataclass()
class InternalStandard():
    """
    Data container to store user input.
    """
    key: str
    conc: Optional[float]

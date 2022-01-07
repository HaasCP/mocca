#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 14:11:19 2022

@author: haascp
"""


def round_to_n(x, n):
    """
    Returns number in a format suitable for data visualization.
    """
    if 1e-3 < x < 1e3:
        return round(x, n)
    else:
        return "{0:.{1}e}".format(x, n)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 09:01:51 2022

@author: haascp
"""


def check_comp_overlap(peak, comp):
    """
    Checks if a given peak overlaps with a given component.
    """
    return comp.left <= peak.left - peak.offset <= comp.right\
        or peak.left - peak.offset <= comp.left <= peak.right - peak.offset


def check_any_compound_overlap(peak, quali_comp_db):
    """
    Checks if a given peak overlaps with any component in the quali_comp_db.
    """
    return any((check_comp_overlap(peak, comp) and 'unknown' not in comp.compound_id
                and 'impurity' not in comp.compound_id) for comp in quali_comp_db)

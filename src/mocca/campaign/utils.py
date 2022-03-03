#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:47:14 2021

@author: haascp
"""
import dill


def save_instance(campaign, path):
    """
    Saves a HplcCampaign object as pkl file.
    """
    with open(path, 'wb') as file:
        dill.dump(campaign, file)


def load_instance(path):
    """
    Loads a pkl file containing an HplcCampaign object.
    """
    with open(path, 'rb') as file:
        dill.load(file)


def check_istd(exp, chrom):
    """
    Checks internal standard condition, ie, if the user gives an istd information
    in the experiment, a corresponding peak has to be found in the chromatogram.
    """
    if exp.istd:
        for istd in exp.istd:
            if not any([peak.compound_id == istd.key for peak in chrom]):
                chrom.bad_data = True
    return chrom

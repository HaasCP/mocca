# flake8: noqa
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 15:17:53 2021

@author: haascp
"""


def read_custom_data(experiment):
    """
    Returns the given custom data without any preprocessing
    """
    if experiment.custom_data is None:
        raise AttributeError("Custom data has to be given if data should be "
                             "processed with custom hplc_system_tag.")
    custom_data = experiment.custom_data
    return custom_data.data, custom_data.time, custom_data.wavelength

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 16:11:33 2021

@author: haascp
"""
import logging

from mocca.campaign.utils import suggest_initialization_runs

def test_suggest_initialization_runs_1():
    n_calib_dict = {'a': 2, 'b': 1}
    max_conc_dict = {'a': 1, 'b': 2}
    calib_dict = suggest_initialization_runs(n_calib_dict, max_conc_dict=max_conc_dict)
    assert calib_dict['a'][0] == 1 and calib_dict['b'][0] == 0
    assert len(calib_dict['a']) == len(calib_dict['b']) and len(calib_dict['a']) == 3
    for i in range(len(calib_dict['a'])):
        non_zero_counter = 0
        for key in calib_dict.keys():
            if calib_dict[key][i] != 0:
                non_zero_counter+=1
        assert non_zero_counter == 1

def test_suggest_initialization_runs_2():
    n_calib_dict = {'a': 2, 'b': 1}
    max_conc_dict = {'a': 1, 'b': 2, 'c':4}
    istd_key = "c"
    n_solvents = 2
    calib_dict = suggest_initialization_runs(n_calib_dict, max_conc_dict=max_conc_dict,
                                             istd_key=istd_key, n_solvents=n_solvents)
    assert len(calib_dict['a']) == len(calib_dict['c']) and len(calib_dict['a']) == 6
    non_zero_counter_list = []
    for i in range(len(calib_dict['a'])):
        non_zero_counter = 0
        for key in calib_dict.keys():
            if calib_dict[key][i] != 0:
                non_zero_counter+=1
        non_zero_counter_list.append(non_zero_counter)
    assert non_zero_counter_list == [0, 0, 1, 2, 2, 2]
    assert all(val == 4 for val in calib_dict['c'][2:])

#logging.warning("{}".format(calib_dict))
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 19:13:39 2021

@author: haascp
"""
import pytest

from mocca.peak.models import BasePeak, PreprocessedPeak

from mocca.campaign.funcs_init import get_max_integral_peak

def test_get_max_integral_peak_1():
    peak_1 = PreprocessedPeak(left=90, right=110, maximum=100, dataset="a", idx=1,
                              pure=False, saturation=False, integral=0.99, offset=0, matches=[])
    peak_2 = PreprocessedPeak(left=190, right=210, maximum=200, dataset="a", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=[])
    peak_3 = PreprocessedPeak(left=290, right=310, maximum=300, dataset="a", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=[])
    peaks = [peak_1, peak_2, peak_3]
    max_peak = get_max_integral_peak(peaks)
    assert max_peak == peak_2

def test_get_max_integral_peak_2():
    peak_1 = PreprocessedPeak(left=90, right=110, maximum=100, dataset="a", idx=1,
                              pure=False, saturation=False, integral=0.99, offset=0, matches=[])
    peak_2 = PreprocessedPeak(left=190, right=210, maximum=200, dataset="a", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=[])
    peak_3 = BasePeak(left=290, right=310, maximum=300)
    peaks = [peak_1, peak_2, peak_3]
    with pytest.raises(AttributeError):
        max_peak = get_max_integral_peak(peaks)
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 08:39:04 2022

@author: HaasCP
"""

from mocca.report.hplc_input import report_hplc_input
from mocca.report.gradient import report_gradients
from mocca.report.chromatograms import report_chroms
from mocca.report.bad_chromatograms import report_bad_chroms
from mocca.report.compound_tracking import report_comp_tracking
from mocca.report.peak_library import report_peak_library
from mocca.report.compound_library import report_comp_library
from mocca.report.calibration_library import report_calib_library
from mocca.report.deconvolution import report_deconvolution


def report(camp, export_path):
    """
    Consolidated report function.
    """
    report_hplc_input(camp.hplc_inputs, export_path)
    report_gradients(camp.hplc_inputs, export_path)
    report_peak_library(camp.peak_db, export_path)
    report_chroms(camp.chroms, camp.settings, export_path)
    report_bad_chroms(camp.chroms, camp.settings, export_path)
    report_comp_tracking(camp.chroms, camp.quali_comp_db, camp.quant_comp_db, export_path)
    report_deconvolution(camp.chroms, export_path)
    report_comp_library(camp.quali_comp_db, export_path)
    report_calib_library(camp.quant_comp_db, export_path)

# -*- coding: utf-8 -*-
"""
Created on Wed May 25 08:39:04 2022

@author: HaasCP
"""

from mocca.report.hplc_input import report_hplc_input
from mocca.report.gradient import report_gradients
from mocca.report.chroms import report_chroms
from mocca.report.bad_chroms import report_bad_chroms
from mocca.report.results import report_runs
from mocca.report.peaks import report_peaks
from mocca.report.quali_comps import report_quali_comps
from mocca.report.quant_comps import report_quant_comps
from mocca.report.parafac import report_parafac


def report(camp, export_path=''):
    """
    Consolidated report function.
    """
    report_hplc_input(camp.hplc_inputs, export_path)
    report_gradients(camp.hplc_inputs, export_path)
    report_peaks(camp.peak_db, export_path)
    report_chroms(camp.chroms, camp.settings, export_path)
    report_bad_chroms(camp.chroms, camp.settings, export_path)
    report_runs(camp.chroms, camp.quali_comp_db, camp.quant_comp_db, export_path)
    report_parafac(camp.chroms, export_path)
    report_quali_comps(camp.quali_comp_db, export_path)
    report_quant_comps(camp.quant_comp_db, export_path)

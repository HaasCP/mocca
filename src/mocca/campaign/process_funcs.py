#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 10:43:13 2021

@author: haascp
"""

from mocca.dad_data.models import CompoundData, GradientData
from mocca.dad_data.process_funcs import pick_peaks

from mocca.chromatogram.preprocessor import preprocess_chromatogram
from mocca.chromatogram.assign import (assign_peaks_compound,
                                       reassign_impurities,
                                       assign_peaks_react)
from mocca.chromatogram.quantify import quantify_peaks

from mocca.campaign.utils import check_istd
from mocca.campaign.experiment_funcs import (get_sorted_compound_experiments,
                                             get_unprocessed_experiments)
from mocca.user_interaction.user_objects import HplcInput


def process_gradients(experiments, settings):
    """
    Reads and processes gradient data for each experiment. Avoids double processing.
    """
    unprocessed_exps = get_unprocessed_experiments(experiments)
    relevant_exps = [exp for exp in unprocessed_exps if exp.gradient]
    for exp in relevant_exps:
        gradient_dataset = next((e.gradient.dataset for e in unprocessed_exps
                                 if hasattr(e.gradient, 'dataset') and
                                 exp.gradient.path == e.gradient.path),
                                None)
        if gradient_dataset is None:
            exp.gradient.dataset = GradientData(settings.hplc_system_tag,
                                                HplcInput(exp.gradient.path, None),
                                                settings.wl_high_pass,
                                                settings.wl_low_pass)
        else:
            exp.gradient.dataset = gradient_dataset


def preprocess_experiment(exp, quali_comp_db, settings):
    """
    Returns a chromatogram object created out of the given experiment.
    The peaks in the chromatogram already have assigned possible matches
    but they are not yet assigned or quantified.
    """
    compound_data = CompoundData(settings.hplc_system_tag, exp,
                                 wl_high_pass=settings.wl_high_pass,
                                 wl_low_pass=settings.wl_low_pass)
    chromatogram = pick_peaks(compound_data, exp, settings.absorbance_threshold,
                              settings.peaks_high_pass,
                              settings.peaks_low_pass)
    chromatogram.experiment = exp
    chromatogram = preprocess_chromatogram(chromatogram, quali_comp_db,
                                           settings.absorbance_threshold,
                                           settings.detector_limit,
                                           settings.spectrum_correl_thresh,
                                           settings.relative_distance_thresh)
    return chromatogram


def process_compound_exp(exp, quali_comp_db, settings):
    """
    Processes one compound experiment.
    """
    chromatogram = preprocess_experiment(exp, quali_comp_db, settings)
    chromatogram = assign_peaks_compound(chromatogram, exp.compound)
    return chromatogram


def process_compound_experiments(experiments, peak_db, quali_comp_db,
                                 quant_comp_db, settings):
    """
    Sorts and processes all compound experiments (experiments from which the
    program learns). Updates both qualitative and quantitative component
    databases.
    """
    exps = get_sorted_compound_experiments(experiments)
    chroms = []
    for exp in exps:
        chrom = process_compound_exp(exp, quali_comp_db, settings)
        if not chrom.bad_data:
            chrom = check_istd(exp, chrom)
        chroms.append(chrom)
        if not chrom.bad_data:
            for peak in chrom:
                if 'impurity' not in peak.compound_id:
                    peak_db.insert_peak(peak)
            quali_comp_db.update(peak_db)

    for chrom in [chrom for chrom in chroms if not chrom.bad_data]:
        chrom = reassign_impurities(chrom, peak_db, quali_comp_db,
                                    settings.spectrum_correl_thresh,
                                    settings.relative_distance_thresh)
        for peak in chrom:
            if peak not in peak_db and peak.idx > 0:
                peak_db.insert_peak(peak)
        quali_comp_db.update(peak_db)

    quant_comp_db.update(peak_db, quali_comp_db)

    return chroms


def process_experiments(experiments, peak_db, quali_comp_db, quant_comp_db,
                        settings):
    """
    Processes all unprocessed experiments (not compound experiments) which should
    be analyzed by the program.
    """
    unprocessed_exps = get_unprocessed_experiments(experiments, quali_comp_db)
    chroms = []
    for exp in unprocessed_exps:
        chrom = preprocess_experiment(exp, quali_comp_db, settings)
        chrom = assign_peaks_react(chrom, peak_db)
        chrom = quantify_peaks(chrom, quant_comp_db, quali_comp_db)
        chrom = check_istd(exp, chrom)
        chroms.append(chrom)
        if not chrom.bad_data:
            for peak in chrom:
                peak_db.insert_peak(peak)
            quali_comp_db.update(peak_db)
    return chroms

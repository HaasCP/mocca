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


def get_gradient_experiment(experiments):
    gradient_experiments = [exp for exp in experiments if exp.gradient]
    warnings = []
    if not gradient_experiments:
        raise ValueError("Gradient run must be provided for the campaign! "
                         "Add blank gradient run data via add_experiment "
                         "function with the gradient attribute equals True.")
        
    elif len(gradient_experiments) > 1:
        warnings.append("Gradient Warning: More than one gradient run were "
                        "detected in the experiments. Per default, the latest "
                        "is taken for this campaign.")
    return gradient_experiments[-1]


def get_gradient(experiments, settings):
    gradient_experiment = get_gradient_experiment(experiments)
    return GradientData(settings.hplc_system_tag, gradient_experiment,
                        settings.wl_high_pass, settings.wl_low_pass)


def get_sorted_compound_experiments(experiments):
    """
    Fitlers experiments for experiments with given compound. Sorts these
    experiments in the order: 1. solvent runs, 2. istd runs, 3. compound runs.
    In these categories, experiments are sorted in the order the user has given.
    """
    compound_exps = [exp for exp in experiments if exp.compound]
    
    solvent_exps = [exp for exp in compound_exps if exp.compound.solvent]

    istd_exps = [exp for exp in compound_exps if exp.compound.istd]
    
    other_exps = [exp for exp in compound_exps if not
                  exp.compound.solvent and not exp.compound.istd]

    conc_exps = [exp for exp in other_exps if exp.compound.conc]
    sorted_conc_exps = sorted(conc_exps, key=lambda exp:
                              (exp.compound.key, -exp.compound.conc))
    non_conc_exps = [exp for exp in other_exps if not exp.compound.conc]

    return solvent_exps + istd_exps + sorted_conc_exps + non_conc_exps


def preprocess_experiment(exp, gradient, quali_comp_db, settings):
    compound_data = CompoundData(settings.hplc_system_tag, exp, gradient,
                                 wl_high_pass=settings.wl_high_pass,
                                 wl_low_pass=settings.wl_low_pass)
    chromatogram = pick_peaks(compound_data, settings.absorbance_threshold, 
                              settings.peaks_high_pass,
                              settings.peaks_low_pass)
    chromatogram = preprocess_chromatogram(chromatogram, exp.istd,
                                           quali_comp_db, 
                                           settings.absorbance_threshold, 
                                           settings.detector_limit, 
                                           settings.spectrum_correl_thresh,
                                           settings.relative_distance_thresh)
    return chromatogram


def process_compound_exp(exp, gradient, quali_comp_db, settings):
    """
    Processes one compound experiment.
    """
    chromatogram = preprocess_experiment(exp, gradient, quali_comp_db, settings)
    chromatogram = assign_peaks_compound(chromatogram, exp.compound)
    return chromatogram


def check_istd(exp, chrom):
    if exp.istd:
        for istd in exp.istd:
            if not any([peak.compound_id == istd.key for peak in chrom]):
                chrom.bad_data = True
    return chrom


def process_compound_experiments(experiments, gradient, peak_db,
                                 quali_comp_db, quant_comp_db, settings):
    exps = get_sorted_compound_experiments(experiments)
    compound_chroms = []
    bad_chroms = []
    for exp in exps:
        chrom = process_compound_exp(exp, gradient,
                                     quali_comp_db, settings)
        chrom = check_istd(exp, chrom)
        if not chrom.bad_data:
            compound_chroms.append(chrom)
            for peak in chrom:
                if not 'impurity' in peak.compound_id:
                    peak_db.insert_peak(peak)
            quali_comp_db.update(peak_db)
        else:
            bad_chroms.append(chrom)

    for chrom in compound_chroms:
        chrom = reassign_impurities(chrom, peak_db, quali_comp_db,
                                    settings.spectrum_correl_thresh,
                                    settings.relative_distance_thresh)
        for peak in chrom:
            if not peak in peak_db and peak.idx > 0:
                peak_db.insert_peak(peak)
        quali_comp_db.update(peak_db)

    quant_comp_db.update(peak_db)

    return compound_chroms, bad_chroms


def process_experiments(experiments, gradient, peak_db, quali_comp_db,
                        quant_comp_db, settings):
    unprocessed_exps = [exp for exp in experiments if not exp.processed]
    chroms = []
    bad_chroms = []
    for exp in unprocessed_exps:
        chrom = preprocess_experiment(exp, gradient, quali_comp_db, settings)
        chrom = assign_peaks_react(chrom, peak_db)
        chrom = quantify_peaks(chrom, quant_comp_db)
        chrom = check_istd(exp, chrom)

        if not chrom.bad_data:
            chroms.append(chrom)
            for peak in chrom:
                    peak_db.insert_peak(peak)
            quali_comp_db.update(peak_db)
        else:
            bad_chroms.append(chrom)    
    return chroms, bad_chroms
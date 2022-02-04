#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 09:11:54 2022

@author: haascp
"""


def get_sorted_compound_experiments(experiments):
    """
    Filters experiments for experiments with given compound. Sorts these
    experiments in the order: 1. solvent runs, 2. istd runs, 3. compound runs
    (sorted reversely by the compound concentration).
    In these categories, experiments are sorted in the order the user has given.
    """
    compound_exps = [exp for exp in experiments if exp.compound]
    
    solvent_exps = [exp for exp in compound_exps if exp.compound.is_solvent]

    istd_exps = [exp for exp in compound_exps if exp.compound.is_istd]
    
    other_exps = [exp for exp in compound_exps if not
                  exp.compound.is_solvent and not exp.compound.is_istd]

    conc_exps = [exp for exp in other_exps if exp.compound.conc]
    sorted_conc_exps = sorted(conc_exps, key=lambda exp:
                              (exp.compound.key, -exp.compound.conc))
    non_conc_exps = [exp for exp in other_exps if not exp.compound.conc]

    return solvent_exps + istd_exps + sorted_conc_exps + non_conc_exps

def get_unprocessed_experiments(experiments, quali_comp_db=None):
    """
    Returns all experiments which have not been processed yet. Checks for
    internal standard condition, ie, that any given istd given in unprocessed
    peaks is already in the qualitative component db so that a corresponding
    peak can be found in the chromatogram.
    """
    unprocessed_exps = [exp for exp in experiments if not exp.processed]
    if quali_comp_db:
        for exp in unprocessed_exps:
            if exp.istd:
                for istd in exp.istd:
                    if istd.key not in quali_comp_db:
                        raise ValueError("Internal standard {} unknown in this campaign. "
                                         "First add the internal standard as pure "
                                         "compound in a separate run!".format(exp.istd.key))
    return unprocessed_exps
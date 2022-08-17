#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 15:44:47 2021

@author: haascp
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def bsl_als_alg(y, lam=1e5, p=0.01, niter=3):
    """
    Baseline correction algorithm: Optimized Python implementation of
    "Asymmetric Least Squares Smoothing" by P. Eilers and H. Boelens in 2005:
    https://stackoverflow.com/questions/29156532/python-baseline-correction-library,
    answer by Rustam Guliev.

    Parameters
    ----------
    y : list
        List of absorbance values for which the baseline should be determined.
    lam : numeric, optional
        Smoothness parameter. 10^2 ≤ λ ≤ 10^9, but exceptions may occur.
        In any case one should vary λ on a grid that is approximately linear
        for log λ. Often visual inspection is sufficient for good values.
        The default is 1e6.
    p : numeric, optional
        Asymmetry parameter. 0.001 ≤ p ≤ 0.1 (for a signal with positive peaks),
        but exceptions may occur. Often visual inspection is sufficient to get
        good parameter values.
        The default is 0.01.
    niter : integer, optional
        To emphasize the basic simplicity of the algorithm, the number of
        iterations has been fixed to 10 (original documentation). In practical
        applications one should check whether the weights show any change;
        if not, convergence has been attained.
        The default is 3.

    Returns
    -------
    z : list
        Simulated baseline of the given absorbance data.
    """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    # Precompute this term since it does not depend on `w`
    D = lam * D.dot(D.transpose())
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w)  # Do not create a new matrix, just update diagonal values
        Z = W + D
        # somehow necessary since matrix can get singular
        # https://stackoverflow.com/questions/48370035/a-is-not-invertible-but-\
        # a-1-00001-is
        Z = Z * 1.00000000000001
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z


def bsl_als(absorbance_array):
    """
    Applies the baseline als algorithm row-wise (for every wavelength) on an
    absorbance array

    Parameters
    ----------
    absorbance_arry : numpy 2D-array
        Absorbance values obtained by an HPLC run (time, wavelength dimension).

    Returns
    -------
    baseline_array : numpy 2D-array
        Baseline absorbance values
    """
    baseline_array = np.apply_along_axis(bsl_als_alg, 1, absorbance_array)
    return baseline_array

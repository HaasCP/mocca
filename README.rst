.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/mocca.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/mocca
    .. image:: https://readthedocs.org/projects/mocca/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://mocca.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/mocca/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/mocca
    .. image:: https://img.shields.io/pypi/v/mocca.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/mocca/
    .. image:: https://img.shields.io/conda/vn/conda-forge/mocca.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/mocca
    .. image:: https://pepy.tech/badge/mocca/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/mocca
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/mocca

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

.. image:: https://github.com/haascp/mocca/blob/master/docs/mocca_icon_w.png?raw=true

|

    MOCCA (Multivariate Online Contextual Chromatographic Analysis) is an open-source Python project to analyze HPLC–DAD raw data.


Automation and digitalization in the field of small molecule synthesis remains a challenge due to the diversity and complexity of chemical reaction processes. Organic chemical reactions are typically developed in Design–Make–Test–Analyze (DMTA) cycles. While much emphasis is given in recent literature on experimental design algorithms (D) and automated reaction execution (M), the analytical side of the cycle, especially the analysis of analytical data (A), remains locked up in the hands of analytical instrument vendors. In this work, we present an open-source Python package called MOCCA (Multivariate Online Contextual Chromatographic Analysis) for the analysis of raw data obtained from high-performance liquid chromatography analysis with diode array detectors (HPLC–DAD), a standard analytical tool to analyze outcomes of small molecule reactions. We investigate the potential of moving HPLC–DAD data analysis into an open environment like Python while overcoming vendor specific proprietary file formats and analysis algorithms. With a focus on automated DMTA workflows, we develop data analysis features, including the checking peak purity and automated deconvolution of overlapping peaks.

Check out corresponding scientific publication:
!!!

Installation
============
#. We recommend creating an isolated conda environment 
   to avoid any problems with your installed Python packages::

    conda create -n mocca python=3.9
    conda activate mocca

#. Install ``mocca``::

    pip install mocca

#. If you want to use ``mocca``'s reporting functionality::

    pip3 install -U datapane

#. If you want to use Allotrope (adf) file format::

    pip install h5py
    pip install git+https://github.com/HDFGroup/h5ld@master

#. If you want to use ``mocca`` using JupyterLab notebooks::

    pip install jupyterlab
    ipython kernel install --user --name=mocca


.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.1.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.

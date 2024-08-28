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

⛔️ DEPRECATED
============

For updated version, please see:

https://github.com/Bayer-Group/MOCCA

|

.. image:: https://github.com/haascp/mocca/blob/master/docs/mocca_icon_w.png?raw=true

|

    MOCCA (Multivariate Online Contextual Chromatographic Analysis) is an open-source Python project to analyze HPLC–DAD raw data.


Automation and digitalization solutions in the field of small molecule synthesis face new challenges for chemical reaction analysis, especially in the field of high-performance liquid chromatography (HPLC). Chromatographic data remains locked in vendors’ hardware and software components limiting their potential in automated workflows and contradicting to FAIR data principles (findability, accessibility, interoperability, reuse), which enable chemometrics and data science applications. In this work, we present an open-source Python project called MOCCA (Multivariate Online Contextual Chromatographic Analysis) for the analysis of open-format HPLC–DAD (photodiode array detector) raw data. MOCCA provides a comprehensive set of data analysis features including a peak deconvolution routine which allows for automated deconvolution of known signals even if overlapped with signals of unexpected impurities or side products. By publishing MOCCA as a Python package, we envision an open-source community project for chromatographic data analysis with the potential of further advancing its scope and capabilities.

Open-source project: https://github.com/HaasCP/mocca

Documentation: https://mocca.readthedocs.io/en/latest/

Corresponding scientific publication (open access): *ACS Central Science*, **2023**, https://doi.org/10.1021/acscentsci.2c01042.

Installation
============
#. We recommend creating an isolated conda environment 
   to avoid any problems with your installed Python packages::

    conda create -n mocca python=3.9
    conda activate mocca

#. Install ``mocca`` and its dependencies::

    pip install mocca

#. If you want to use ``mocca``'s reporting functionality::

    pip3 install -U datapane==0.14

#. If you want to use Allotrope (adf) file format::

    pip install h5py
    pip install git+https://github.com/HDFGroup/h5ld@master

#. If you want to use ``mocca`` using JupyterLab notebooks::

    pip install jupyterlab
    ipython kernel install --user --name=mocca


Getting started
===============
MOCCA is currently best used via JupyterLab notebooks. The notebooks folder of the GitHub repository contains a tutorial notebook with corresponding HPLC–DAD test data for the first steps.

Additionally, a full test data set from the scientific publication is added (cyanation of aryl halides via well plate screening). The corresponding notebook contains full data analysis details from the raw data level until the presented visualizations in the manuscript (Fig. 7e) and SI (Fig. S17).


How to cite
===========
**Peer-reviewed, open access:**

Haas, C. P., Lübbesmeyer, M., Jin, E. H., McDonald, M. A., Koscher, B. A., Guimond, N., Di Rocco, L., Kayser, H., Leweke, S., Niedenführ, S., Nicholls, R., Greeves, E., Barber, D. M., Hillenbrand, J., Volpin, G., and Jensen, K. F. Open-Source Chromatographic Data Analysis for Reaction Optimization and Screening. *ACS Cent.Sci.* **2023**. https://doi.org/10.1021/acscentsci.2c01042.

**Preprint:**

Haas, C. P., Lübbesmeyer, M., Jin, E. H., McDonald, M. A., Koscher, B. A., Guimond, N., Di Rocco, L., Kayser, H., Leweke, S., Niedenführ, S., Nicholls, R., Greeves, E., Barber, D. M., Hillenbrand, J., Volpin, G., and Jensen, K. F. Open-Source Chromatographic Data Analysis for Reaction Optimization and Screening. *ChemRxiv* **2022**. https://doi.org/10.26434/chemrxiv-2022-0pv2d.


.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.1.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.

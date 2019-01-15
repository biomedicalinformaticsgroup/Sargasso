[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597578.svg)](https://doi.org/10.5281/zenodo.597578)
[![Build Status](https://travis-ci.org/statbio/Sargasso.svg?branch=master)](https://travis-ci.org/statbio/Sargasso)

Sargasso
========

*Sargasso* is a Python tool to disambiguate mixed-species high-throughput sequencing reads according to their species of origin. Given a set of samples containing sequencing data from multiple species, mapped, disambiguated reads are written to per-sample and species-specific output BAM files.

The latest *Sargasso* documentation can be found [here](http://statbio.github.io/Sargasso/).

Installation
============

Note: as *Sargasso* has a number of dependencies on other Python packages, it is **strongly** recommended to install in an isolated environment using the [virtualenv](http://virtualenv.readthedocs.org/en/latest/index.html>) tool. The [virtualenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/install.html>) tool makes managing multiple virtual environments easier.

After setting up ``virtualenv`` and ``virtualenvwrapper``, create and work in a virtual environment for *Sargasso* using the ``virtualenvwrapper`` tool:

```
mkproject sargasso
```

Then install the *Sargasso* package and its Python package dependencies into the virtual environment by running:

```
pip install git+https://github.com/statbio/sargasso.git
```

Note that Sargasso v2.0 should work with Python versions >= 2.6 (including Python 3). Versions before v1.2.2 will only work with Python 2.

Citation
========

If you make use of *Sargasso* please cite [our protocol paper](https://www.nature.com/articles/s41596-018-0029-2):

* Qiu *et al.*, "Mixed-species RNA-seq for elucidation of non-cell-autonomous control of gene transcription", *Nature Protocols* **13**, 2176–2199 (2018).

Changelog
=========

* 2.0 (16/01/2019): Sargasso now separates reads derived from DNA-based sequencing technologies (for example, ChIP-seq and ATAC-seq), in addition to RNA-seq reads.
* 1.2.2 (11/10/2018): Bugfix release for compatibility with Python 3.
* 1.2.1 (02/10/2018): Bugfix release for incompatibilities between Mac OS and Linux.
* 1.2 (16/02/2018): 
    - Improvements to species read assignment logic gives better precision and recall.
    - Added --delete-intermediate option to delete intermediate BAM files.
    - Added --star-executable option to allow different versions of STAR to be used.
    - Added --sambamba-sort-tmp-dir option to specify a different temporary directory for 'sambamba sort'.
* 1.1.2 (14/12/2017): Minor improvements to interpretability of results.
* 1.1.1 (02/03/2017): Add "permissive" filtering strategy.
* 1.1 (26/01/2017): Filtering of RNA-seq data from more than two species.
* 1.0 (16/12/2016): First full release

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.206619.svg)](https://doi.org/10.5281/zenodo.206619)
[![Build Status](https://travis-ci.org/statbio/Sargasso.svg?branch=master)](https://travis-ci.org/statbio/Sargasso)

Sargasso
========

*Sargasso* is a Python tool to disambiguate mixed-species RNA-seq reads according to their species of origin. Given a set of RNA-seq samples containing RNA-seq data originating from two different species, mapped, disambiguated reads are written to per-sample and -species specific output BAM files.

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

Changelog
=========

* 1.1 (26/01/2017): Filtering of RNA-seq data from more than two species.
* 1.0 (16/12/2016): First full release

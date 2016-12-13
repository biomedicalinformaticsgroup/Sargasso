Sargasso
========

*Sargasso* is a Python tool to disambiguate mixed-species RNA-seq reads according to their species of origin. Given a set of RNA-seq samples containing RNA-seq data originating from two different species, mapped, disambiguated reads are written to per-sample and -species specific output BAM files.

The latest *Sargasso* documentation can be found **docs placeholder**.

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

* 1.0 (xx/xx/16): First full release

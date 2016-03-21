Species Separator
=================

*species-separator* is a Python tool to disambiguate mixed-species RNA-seq reads according to their species of origin. Given a set of RNA-seq samples containing RNA-seq data originating from two different species, mapped, disambiguated reads are written to per-sample and -species specific output BAM files.

The latest *species-separator* documentation can be found **docs placeholder**.

Installation
============

Note: as *species-separator* has a number of dependencies on other Python packages, it is **strongly** recommended to install in an isolated environment using the [virtualenv](http://virtualenv.readthedocs.org/en/latest/index.html>) tool. The [virtualenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/install.html>) tool makes managing multiple virtual environments easier.

After setting up ``virtualenv`` and ``virtualenvwrapper``, create and work in a virtual environment for *species-separator* using the ``virtualenvwrapper`` tool:

```
mkproject species_separator
```

Install the *species-separator* package and its Python package dependencies into the virtual environment by running:

```
pip install git+https://github.com/lweasel/species-separator.git
```

Changelog
=========

* 1.0 (xx/xx/16): First full release

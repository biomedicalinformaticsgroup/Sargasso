Installation
============

As *Sargasso* has a number of dependencies on other Python packages, it is **strongly** recommended to install in an isolated environment using the [virtualenv](http://virtualenv.readthedocs.org/en/latest/index.html) tool. The [virtualenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/install.html) tool makes managing multiple virtual environments easier.

Create and work in a virtual environment for *Sargasso* using the ``virtualenvwrapper`` tool::

    mkproject sargasso

Then install the *Sargasso* package and its Python package dependencies into the virtual environment by running::

    pip install git+https://github.com/statbio/Sargasso.git

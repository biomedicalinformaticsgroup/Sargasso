Installation
============

.. attention:: As *species-separator* has a number of dependencies on other Python packages, it is **strongly** recommended to install in an isolated environment using the `virtualenv <http://virtualenv.readthedocs.org/en/latest/index.html>`_ tool. The `virtualenvwrapper <http://virtualenvwrapper.readthedocs.org/en/latest/install.html>`_ tool makes managing multiple virtual environments easier.

Create and work in a virtual environment for *species-separator* using the ``virtualenvwrapper`` tool::

    mkproject species_separator

Install the *species_separator* package and its Python package dependencies into the virtual environment by running::

    pip install git+https://github.com/statbio/species-separator.git

.. Tests of the species separation pipeline can be run by executing::

..    (cd pipeline_test; ./run_test.sh)

.. These tests should take a couple of minutes to run; no output, with exit code 0, indicates successful test completion.

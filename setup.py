from setuptools import setup

import species_separator

setup(
    name="species_separator",
    version=species_separator.__version__,
    url='https://github.com/statbio/Sargasso',
    license='MIT License',
    author='Owen Dando',
    author_email='owen.dando@ed.ac.uk',
    packages=['species_separator'],
    install_requires=[
        'docopt==0.6.2',
        'pysam==0.8.2.1',
        'schema==0.3.1',
    ],
    scripts=[
        'bin/build_star_index',
        'bin/collate_raw_reads',
        'bin/filter_control',
        'bin/filter_reads',
        'bin/filter_sample_reads',
        'bin/map_reads',
        'bin/sort_reads',
        'bin/species_separator',
    ]
)

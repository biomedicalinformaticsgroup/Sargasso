from setuptools import setup, find_packages

import sargasso

setup(
    name="sargasso",
    version=sargasso.__version__,
    url='https://github.com/statbio/Sargasso',
    license='MIT License',
    author='Owen Dando',
    author_email='owen.dando@ed.ac.uk',
    packages=find_packages(),
    install_requires=[
        'docopt',
        'pysam',
        'schema',
        'pytest',
    ],
    scripts=[
        'bin/build_star_index',
        'bin/build_bowtie2_index',
        'bin/collate_raw_reads',
        'bin/filter_control',
        'bin/filter_reads',
        'bin/filter_sample_reads',
        'bin/map_reads_rnaseq',
        'bin/map_reads_dnaseq',
        'bin/sort_reads',
        'bin/species_separator',
    ]
)

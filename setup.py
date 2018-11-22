from setuptools import setup, find_packages

import sargasso

setup(
    name="sargasso",
    version=sargasso.__version__,
    url='https://github.com/statbio/Sargasso',
    license='MIT License',
    author='Owen Dando',
    author_email='owen.dando@ed.ac.uk',
    # packages=['sargasso'],
    packages=find_packages(),
    install_requires=[
        'docopt==0.6.2',
        'pysam==0.8.2.1',
        'schema==0.3.1',
        'pytest==3.9.1',
    ],
    scripts=[
        'bin/build_star_index',
        'bin/build_bowtie2_index',
        'bin/collate_raw_reads',
        'bin/filter_control',
        'bin/filter_reads',
        'bin/filter_sample_reads',
        'bin/map_reads_rnaseq',
        'bin/map_reads_chipseq',
        'bin/sort_reads',
        'bin/species_separator',
    ]
)

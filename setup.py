from setuptools import setup

import species_separator

setup(
    name="species_separator",
    version=species_separator.__version__,
    url='https://github.com/lweasel/species-separator',
    license='MIT License',
    author='TODO',
    author_email='TODO',
    packages=['species_separator'],
    install_requires=[
        'docopt==0.6.2',
        'pysam==0.8.2.1',
        'schema==0.3.1',
    ],
    scripts=[
        'bin/build_star_index',
        'bin/collate_raw_reads_pe',
        'bin/collate_raw_reads_se',
        'bin/filter_control',
        'bin/filter_reads',
        'bin/filter_sample_reads',
        'bin/map_reads_pe',
        'bin/map_reads_se',
        'bin/sort_reads',
        'bin/species_separator',
    ]
)

from setuptools import setup

import species_separator

setup(
    name="species_separator",
    version=species_separator.__version__,
    url='https://github.com/lweasel/species-separator',
    license='MIT License',
    author='Owen Dando',
    author_email='owen.dando@ed.ac.uk',
    packages=['species_separator'],
    install_requires=[
        'Babel==2.2.0',
        'Jinja2==2.8',
        'MarkupSafe==0.23',
        'Pygments==2.1.3',
        'Sphinx==1.4',
        'alabaster==0.7.7',
        'docopt==0.6.2',
        'docutils==0.12',
        'imagesize==0.7.0',
        'pysam==0.8.2.1',
        'pytz==2016.3',
        'schema==0.3.1',
        'six=1.10.0',
        'snowballstemmer==1.2.1',
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

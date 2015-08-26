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
        'schema==0.3.1',
    ],
    scripts=[
        'bin/filter_mapped_hits',
    ]
)

#!/usr/bin/env python

"""Usage:
    species_separator rnaseq [options]
        <samples-file> <output-dir>
        (<species> <species-star-info>)
        (<species> <species-star-info>)
        ...
    species_separator chipseq [options]
        <samples-file> <output-dir>
        (<species> <species-star-info>)
        (<species> <species-star-info>)
        ...

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<samples-file>
    TSV file giving paths (relative to <reads-base-dir>) of raw RNA-seq read
    data files for each sample.
<output-dir>
    Output directory into which Makefile will be written, and in which species
    separation will be performed.
<species>
    Name of species.
<species-star-info>
    Either a STAR index directory for the species, or a comma-separated list of
    (i) a GTF annotation file and (ii) a directory containing genome FASTA
    files for the species.
<species-genome-fasta>
    Directory containing genome FASTA files for species.
<species-index>
    STAR index directory for species.
"""

import docopt
import os
import os.path
import schema

# from . import options as opt
from .__init__ import __version__
from separators import SeparatorManager
from separators import Separator
from commandline_parser import CommandlineParser
from parameter_validator import ParameterValidator
from options import Options


def separate_species(args):
    """
    Write a Makefile to perform species separation and optionally execute it.

    args: list of command-line arguments
    """
    # Read in command-line options
    # docstring = opt.substitute_common_options_into_usage(__doc__)
    # options = docopt.docopt(docstring, argv=args,
    #                         version="species_separator v" + __version__)

    #https://github.com/docopt/docopt/blob/master/examples/git/git.py
    ops = CommandlineParser().parse_extra(args,Separator.DOC,options_first=True)

    # todo this could be refactor into validator somehow?
    data_type=ops[Options.DATA_TYPE]
    try:
        ParameterValidator.validate_datatype(data_type)
    except schema.SchemaError as exc:
        exit("Exiting: " + exc.code)

    separator = SeparatorManager().get(data_type)
    separator.run(args)

    # todo remove debug code
    print("--------------------------")
    print("make file content:")
    print("--------------------------")

    with open('/home/xinhe/Projects/Sargasso/results/chipseqtest/Makefile', 'r') as fin:
        print fin.read()
    print("--------------------------")

    print('done')
    exit(1)

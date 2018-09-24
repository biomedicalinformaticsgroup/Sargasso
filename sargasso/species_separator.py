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

from . import options as opt
from .__init__ import __version__
from . import separators as ss


def separate_species(args):
    """
    Write a Makefile to perform species separation and optionally execute it.

    args: list of command-line arguments
    """
    # Read in command-line options
    # docstring = opt.substitute_common_options_into_usage(__doc__)
    # options = docopt.docopt(docstring, argv=args,
    #                         version="species_separator v" + __version__)


    separator = ss.SeparatorManager(args[0],args).get_separator()
    separator.run()

    # todo remove debug code
    print("--------------------------")
    print("make file content:")
    print("--------------------------")

    with open('/home/xinhe/Projects/Sargasso/results/anothertest4/Makefile', 'r') as fin:
        print fin.read()
    print("--------------------------")

    print('done')
    exit(1)
    # Validate command-line options
    sample_info = _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Create output directory
    os.mkdir(options[OUTPUT_DIR])

    # Write Makefile to output directory
    _write_makefile(logger, options, sample_info)

    # Write Execution Record to output directory
    _write_execution_record(options)

    # If specified, execute the Makefile with nohup
    if options[RUN_SEPARATION]:
        _run_species_separation(logger, options)

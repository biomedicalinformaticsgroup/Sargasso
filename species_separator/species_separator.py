#!/usr/bin/env python

"""Usage:
    species_separator [--log-level=<log-level>] [--s1-gtf=<species-one-gtf-file>] [--s2-gtf=<species-two-gtf-file>] [--s1-index=<species-one-star-index>] [--s2-index=<species-two-star-index>] [--run-separation] <species-one> <species-two> <samples-file> <output-dir>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<species-one>                           Name of first species.
<species-two>                           Name of second species.
<samples-file>                          TSV file giving raw RNA-seq data files for each sample.
<output-dir>                            Output directory in which species separation will be performed.
--s1-gtf=<species-one-gtf-file>         GTF annotation file for first species.
--s2-gtf=<species-two-gtf-file>         GTF annotation file for second species.
--s1-index=<species-one-star-index>     STAR index directory for first species.
--s2-index=<species-two-star-index>     STAR index directory for first species.
--run-separation                        If specified, species separation will be run; otherwise scripts to perform separation will be created but not run.

TODO: what does this script do...
"""

import docopt

from . import options as opt
from .__init__ import __version__


def _validate_command_line_options(options):
    # validate log level
    # Check samples TSV file exists
    # Check output directory does not already exist
    # Check GTF files for each species exist if specified
    # Check STAR index directories for each species exist if specified
    pass


def _separate_species(logger, options):
    # Create output directory
    # Write Makefile to output directory, which, when executed, will perform species separation
    # If --run-separation has been specified, execute this Makefile (with nohup)
    pass


def separate_species(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="species_separator v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    _separate_species(logger, options)

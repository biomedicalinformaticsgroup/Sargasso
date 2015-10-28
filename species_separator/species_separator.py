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
import os

from . import file_writer as fw
from . import options as opt
from .__init__ import __version__

OUTPUT_DIR = "<output-dir>"
RUN_SEPARATION = "--run-separation"


def _validate_command_line_options(options):
    # TODO: validate log level
    # TODO: Check samples TSV file exists
    # TODO: Check output directory does not already exist
    # TODO: Check GTF files for each species exist if specified
    # TODO: Check STAR index directories for each species exist if specified
    pass


def _write_variable_definitions(logger, writer, options):
    # TODO: Write variable definitions
    writer.add_line("NUM_THREADS=8")


def _write_target_variable_definitions(logger, writer, options):
    # TODO: Write target variable definitions
    pass


def _write_phony_targets(logger, writer, options):
    # TODO: Write .PHONY targets
    pass


def _write_all_target(logger, writer, options):
    # TODO: Write "all" target
    pass


def _write_filtered_reads_target(logger, writer, options):
    # TODO: Write "filtered_reads" target
    pass


def _write_sorted_reads_target(logger, writer, options):
    # TODO: Write "sorted_reads" target
    pass


def _write_mapped_reads_target(logger, writer, options):
    # TODO: Write "mapped_reads" target
    pass


def _write_masked_reads_target(logger, writer, options):
    # TODO: Write "masked_reads" target
    pass


def _write_collate_raw_reads_target(logger, writer, options):
    # TODO: Write "collate_raw_reads" target
    pass


def _write_mask_star_index_targets(logger, writer, options):
    # TODO: Write targets for STAR indices for mask sequences
    pass


def _write_main_star_index_targets(logger, writer, options):
    # TODO: Write targets to link to main STAR indices (if we're using
    # pre-built indices) or build main STAR indices (if not)
    pass


def _write_clean_target(logger, writer, options):
    # TODO: Write "clean" target
    pass


def _write_makefile(logger, options):
    # TODO: Write Makefile to output directory, which, when executed, will perform species separation
    with fw.writing_to_file(fw.MakefileWriter, options[OUTPUT_DIR], "Makefile") as writer:
        _write_variable_definitions(logger, writer, options)
        _write_target_variable_definitions(logger, writer, options)
        _write_phony_targets(logger, writer, options)
        _write_all_target(logger, writer, options)
        _write_filtered_reads_target(logger, writer, options)
        _write_sorted_reads_target(logger, writer, options)
        _write_mapped_reads_target(logger, writer, options)
        _write_masked_reads_target(logger, writer, options)
        _write_collate_raw_reads_target(logger, writer, options)
        _write_mask_star_index_targets(logger, writer, options)
        _write_main_star_index_targets(logger, writer, options)
        _write_clean_target(logger, writer, options)


def _run_species_separation(logger, options):
    # TODO: Execute Makefile with nohup
    pass


def _separate_species(logger, options):
    os.mkdir(options[OUTPUT_DIR])

    _write_makefile(logger, options)

    if options[RUN_SEPARATION]:
        _run_species_separation(logger, options)


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

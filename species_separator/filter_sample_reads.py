#!/usr/bin/env python

"""Usage:
    filter_mapped_hits [--log-level=<log-level>] <species-one> <species-one-input-bam> <species-one-output-bam> <species-two> <species-two-input-bam> <species-two-output-bam> <mismatch-threshold> <minmatch-threshold> <multimap-threshold>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<species-one>               Name of first species.
<species-one-input-bam>     BAM file containing read hits against first species' genome.
<species-one-output-bam>    BAM file to which reads assigned to first species after filtering will be written.
<species-two>               Name of first species.
<species-two-input-bam>     BAM file containing read hits against first species' genome.
<species-two-output-bam>    BAM file to which reads assigned to first species after filtering will be written.
<mismatch-threshold>	    Maximum number of mismatches allowed during filtering
<minmatch-threshold>        Minimum number of read bases that must be perfectly matched
<multimap-threshold>        Maximum number of multiple mappings allowed during filtering

TODO: what does this script do...
"""

import docopt
import os.path
import schema

from . import filterer
from . import hits_checker
from . import options as opt
from .__init__ import __version__

SPECIES_ONE = "<species-one>"
SPECIES_ONE_INPUT_BAM = "<species-one-input-bam>"
SPECIES_ONE_OUTPUT_BAM = "<species-one-output-bam>"
SPECIES_TWO = "<species-two>"
SPECIES_TWO_INPUT_BAM = "<species-two-input-bam>"
SPECIES_TWO_OUTPUT_BAM = "<species-two-output-bam>"
MISMATCH_THRESHOLD = "<mismatch-threshold>"
MINMATCH_THRESHOLD = "<minmatch-threshold>"
MULTIMAP_THRESHOLD = "<multimap-threshold>"


def validate_threshold_options(
        options, mismatch_opt_name, minmatch_opt_name, multimap_opt_name):

    options[mismatch_opt_name] = opt.validate_int_option(
        options[mismatch_opt_name],
        "Maximum number of mismatches must be a non-negative integer.",
        0, True)
    options[minmatch_opt_name] = opt.validate_int_option(
        options[minmatch_opt_name],
        "Maximum number of not perfect matches must be a non-negative " +
        "integer.", 0, True)
    options[multimap_opt_name] = opt.validate_int_option(
        options[multimap_opt_name],
        "Maximum number of multiple mappings must be a positive integer.",
        1, True)


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)

        opt.validate_file_option(
            options[SPECIES_ONE_INPUT_BAM],
            "Could not find input BAM file for species 1")
        opt.validate_file_option(
            options[SPECIES_TWO_INPUT_BAM],
            "Could not find input BAM file for species 2")

        validate_threshold_options(options, MISMATCH_THRESHOLD,
                                   MINMATCH_THRESHOLD, MULTIMAP_THRESHOLD)
    except schema.SchemaError as exc:
        exit(exc.code)


def _filter_sample_reads(logger, options):
    logger.info("Starting species separation.")

    h_check = hits_checker.HitsChecker(
        options[MISMATCH_THRESHOLD], options[MINMATCH_THRESHOLD],
        options[MULTIMAP_THRESHOLD], logger)

    s1_filterer = filterer.Filterer(
        1, options[SPECIES_ONE_INPUT_BAM], options[SPECIES_ONE_OUTPUT_BAM],
        h_check, logger)
    s2_filterer = filterer.Filterer(
        2, options[SPECIES_TWO_INPUT_BAM], options[SPECIES_TWO_OUTPUT_BAM],
        h_check, logger)

    s1_read_name = None
    s2_read_name = None

    while True:
        # Attempt to read all hits for the next read in each species' input BAM
        # file
        try:
            s1_read_name = s1_filterer.get_next_read_name()
        except StopIteration:
            # If no more reads are available for species 1, all remaining
            # reads in the input file for species 2 can be written to the
            # output file for species 2, or discarded as ambiguous
            s2_filterer.check_and_write_hits_for_remaining_reads()
            break

        try:
            s2_read_name = s2_filterer.get_next_read_name()
        except StopIteration:
            # If no more reads are available for species 2, all remaining
            # reads in the input file for species 1 can be written to the
            # output file for species 1, or discarded as ambiguous
            s1_filterer.check_and_write_hits_for_remaining_reads()
            break

        if s1_read_name == s2_read_name:
            # Compare the read names; if they are the same, compare the hits
            # for each species to determine which species to assign the read
            # to.
            s1_filterer.compare_and_write_hits(s2_filterer)
        elif s1_read_name < s2_read_name:
            # if species 1 read name < species 2 read name, there are no hits
            # for the species 1 read in species 2; write the hits for this read
            # to the output file for species 1, or discard as ambiguous
            s1_filterer.check_and_write_hits_for_read()
        else:
            # if species 2 read name < species 1 read name, there are no hits
            # for the species 2 read in species 1; write the hits for this read
            # to the output file for species 2, or discard as ambiguous
            s2_filterer.check_and_write_hits_for_read()

    s1_filterer.log_stats()
    s2_filterer.log_stats()

    write_stats(s1_filterer.stats, s2_filterer.stats, options[SPECIES_ONE_OUTPUT_BAM])


# write filter stats to table in file
def write_stats(s1, s2, outBam):
    stats = [
        s1.hits_written, s1.reads_written,
        s1.hits_rejected, s1.reads_rejected,
        s1.hits_ambiguous, s1.reads_ambiguous,
        s2.hits_written, s2.reads_written,
        s2.hits_rejected, s2.reads_rejected,
        s2.hits_ambiguous, s2.reads_ambiguous
    ]

    out_file = os.path.join(
        os.path.dirname(outBam), "filtering_result_summary.txt")
    with open(out_file, 'a') as outf:
        outf.write("\t".join([str(s) for s in stats]) + "\n")


def filter_sample_reads(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="filter_sample_reads v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    _filter_sample_reads(logger, options)

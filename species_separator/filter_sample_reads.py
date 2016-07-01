
#!/usr/bin/env python

"""Usage:
    filter_sample_reads
        [--log-level=<log-level>] [--reject-multimaps]
        <species-one> <species-one-input-bam> <species-one-output-bam>
        <species-two> <species-two-input-bam> <species-two-output-bam>
        <mismatch-threshold> <minmatch-threshold> <multimap-threshold> <overhang-threshold>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<species-one>
    Name of first species.
<species-one-input-bam>
    BAM file containing reads mapped against first species' genome.
<species-one-output-bam>
    BAM file to which read mappings assigned to first species after filtering
    will be written.
<species-two>
    Name of second species.
<species-two-input-bam>
    BAM file containing reads mapped against second species' genome.
<species-two-output-bam>
    BAM file to which read mappings assigned to second species after filtering
    will be written.
--mismatch-threshold=<mismatch-threshold>
    Maximum percentage of read bases allowed to be mismatches against the
    genome during filtering.
--minmatch-threshold=<minmatch-threshold>
    Maximum percentage of read length allowed to not be mapped during
    filtering.
<multimap-threshold>
    Maximum number of multiple mappings allowed during filtering.
--reject-multimaps
    If set, any read which multimaps to *either* species' genome will be
    rejected and not be assigned to either species.
<overhang-threshold>
    The minimum number of bases that are allowed on
    either side of an exon boundary for a read mapping to be accepted

filter_sample_reads takes two BAM files as input, the results of mapping a set
of mixed species RNA-seq reads against the two species' genomes, and
determines, if possible, from which species each read originates. Disambiguated
read mappings are written to two species-specific output BAM files.

In normal operation, the user should not need to execute this script by hand
themselves.

Note: the input BAM files MUST be sorted in read name order. Failure to
ensure input BAM files are correctly sorted will result in erroneous output.
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
REJECT_MULTIMAPS = "--reject-multimaps"
OVERHANG_THRESHOLD = "<overhang-threshold>"


def validate_threshold_options(
        options, mismatch_opt_name, minmatch_opt_name, multimap_opt_name,
        overhang_opt_name):

    options[mismatch_opt_name] = opt.validate_float_option(
        options[mismatch_opt_name],
        "Maximum percentage of mismatches must be a float between 0 and 100",
        0, 100, True)
    options[minmatch_opt_name] = opt.validate_float_option(
        options[minmatch_opt_name],
        "Maximum percentage of read length which does not match must be a " +
        "float between 0 and 100", 0, 100, True)
    options[multimap_opt_name] = opt.validate_int_option(
        options[multimap_opt_name],
        "Maximum number of multiple mappings must be a positive integer",
        1, True)
    options[overhang_opt_name] = opt.validate_int_option(
        options[overhang_opt_name],
        "Minimum overhang threshold must be an integer 0 or greater",
        0, True)


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
                                   MINMATCH_THRESHOLD, MULTIMAP_THRESHOLD, OVERHANG_THRESHOLD)

    except schema.SchemaError as exc:
        exit(exc.code)


def _filter_sample_reads(logger, options):
    logger.info("Starting species separation.")

    h_check = hits_checker.HitsChecker(
        options[MISMATCH_THRESHOLD], options[MINMATCH_THRESHOLD],
        options[MULTIMAP_THRESHOLD], options[REJECT_MULTIMAPS], options[OVERHANG_THRESHOLD], logger)

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

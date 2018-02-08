#!/usr/bin/env python

"""Usage:
    filter_sample_reads
        [--log-level=<log-level>] [--reject-multimaps]
        <mismatch-threshold> <minmatch-threshold> <multimap-threshold>
        (<species> <species-input-bam> <species-output-bam>)
        (<species> <species-input-bam> <species-output-bam>) ...

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<species>
    Name of species.
<species-input-bam>
    BAM file containing reads mapped against species' genome.
<species-output-bam>
    BAM file to which read mappings assigned to species after filtering
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

filter_sample_reads takes a set of BAM files as input, the results of mapping a set
of mixed species RNA-seq reads against a number of species' genomes, and
determines, if possible, from which species each read originates. Disambiguated
read mappings are written to species-specific output BAM files.

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

SPECIES = "<species>"
SPECIES_INPUT_BAM = "<species-input-bam>"
SPECIES_OUTPUT_BAM = "<species-output-bam>"
MISMATCH_THRESHOLD = "<mismatch-threshold>"
MINMATCH_THRESHOLD = "<minmatch-threshold>"
MULTIMAP_THRESHOLD = "<multimap-threshold>"
REJECT_MULTIMAPS = "--reject-multimaps"


def validate_threshold_options(
        options, mismatch_opt_name, minmatch_opt_name, multimap_opt_name):

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


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)

        for index, species in enumerate(options[SPECIES]):
            opt.validate_file_option(
                options[SPECIES_INPUT_BAM][index],
                "Could not find input BAM file for species {i}".format(i=index))

        validate_threshold_options(
            options, MISMATCH_THRESHOLD, MINMATCH_THRESHOLD, MULTIMAP_THRESHOLD)

    except schema.SchemaError as exc:
        exit(exc.code)


# write filter stats to table in file
def _write_stats(filterers, out_bam):

    stats = []

    for filt in filterers:
        fstats = filt.stats

        stats += [fstats.hits_written, fstats.reads_written,
                  fstats.hits_rejected, fstats.reads_rejected,
                  fstats.hits_ambiguous, fstats.reads_ambiguous]

    out_file = os.path.join(
        os.path.dirname(out_bam), "filtering_result_summary.txt")
    with open(out_file, 'a') as outf:
        outf.write("\t".join([str(s) for s in stats]) + "\n")


def _get_next_read_name(filterer):
    read_name = None

    try:
        read_name = filterer.get_next_read_name()
    except StopIteration:
        pass

    return read_name


def _filter_sample_reads(logger, options):
    logger.info("Starting species separation.")

    h_check = hits_checker.HitsChecker(
        options[MISMATCH_THRESHOLD], options[MINMATCH_THRESHOLD],
        options[MULTIMAP_THRESHOLD], options[REJECT_MULTIMAPS], logger)

    filterers = [filterer.Filterer(i + 1, options[SPECIES_INPUT_BAM][i],
                                   options[SPECIES_OUTPUT_BAM][i], logger)
                 for i, species in enumerate(options[SPECIES])]

    all_filterers = filterers

    while True:
        # Retain only filterers which have hits for the remaining reads
        filterers = [f for f in filterers if _get_next_read_name(f) is not None]

        # If no filterers remain, we're done
        if len(filterers) == 0:
            break

        # If only one filterer remains, all remaining reads in the input file
        # for that species can be written to the output file for that species
        # (or discarded as ambiguous, if necessary).
        if len(filterers) == 1:
            h_check.check_and_write_hits_for_remaining_reads(filterers[0])
            break

        # Compare read names for the filterers, and find the set that have the
        # "lowest" name
        competing_filterers = [filterers[0]]
        min_read_name = _get_next_read_name(filterers[0])

        for cfilt in filterers[1:]:
            read_name = _get_next_read_name(cfilt)
            if read_name == min_read_name:
                competing_filterers.append(cfilt)
            elif read_name < min_read_name:
                competing_filterers = [cfilt]
                min_read_name = read_name

        # If there's only one filterer for this read, write hits for that read
        # to the output file for that species (or discard as ambiguous)
        if len(competing_filterers) == 1:
            h_check.check_and_write_hits_for_read(competing_filterers[0])
            continue

        # Otherwise compare the hits for each species to determine which
        # species to assign the read to.
        h_check.compare_and_write_hits(competing_filterers)

    for filt in all_filterers:
        filt.log_stats()

    _write_stats(all_filterers, options[SPECIES_OUTPUT_BAM][0])


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

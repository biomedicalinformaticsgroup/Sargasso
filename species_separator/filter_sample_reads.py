#!/usr/bin/env python

"""Usage:
    filter_mapped_hits [--log-level=<log-level>] <species-one> <species-one-input-bam> <species-one-output-bam> <species-two> <species-two-input-bam> <species-two-output-bam>

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

TODO: what does this script do...
"""

import docopt
import schema

from . import options as opt
from . import samutils
from .__init__ import __version__

SPECIES_ONE = "<species-one>"
SPECIES_ONE_INPUT_BAM = "<species-one-input-bam>"
SPECIES_ONE_OUTPUT_BAM = "<species-one-output-bam>"
SPECIES_TWO = "<species-two>"
SPECIES_TWO_INPUT_BAM = "<species-two-input-bam>"
SPECIES_TWO_OUTPUT_BAM = "<species-two-output-bam>"


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)

        opt.validate_file_option(
            options[SPECIES_ONE_INPUT_BAM],
            "Could not find input BAM file for species 1")
        opt.validate_file_option(
            options[SPECIES_TWO_INPUT_BAM],
            "Could not find input BAM file for species 2")
    except schema.SchemaError as exc:
        exit(exc.code)


def _sanity_check_hits(hits):
    # TODO: check all hits have the same alignment score
    # TODO: check all hits have the same multimap value
    return hits


def _hits_generator(bam_file):
    all_hits = samutils.open_samfile_for_read(bam_file)
    last_hit_name = None
    current_hits = []

    for hit in samutils.all_reads(all_hits):
        if last_hit_name is None:
            last_hit_name = hit.query_name

        if hit.query_name < last_hit_name:
            # TODO: throw an exception if hits are out of read name order, and
            # handle this gracefully
            pass
        elif hit.query_name == last_hit_name:
            current_hits.append(hit)
        else:
            yield _sanity_check_hits(current_hits)
            last_hit_name = hit.query_name
            current_hits = [hit]

    yield current_hits


def _check_hits_for_read(hits_for_read):
    # TODO: check that the hits for a read are - in themselves - satisfactory
    # to be unambiguously assigned to one species or another. For example, in
    # the ultra-conservative strategy, this means:
    # - no multi-mapping
    # - no mismatches
    # - CIGAR string satisfactory
    return True


def _write_hits(hits_for_read, output_bam):
    if _check_hits_for_read(hits_for_read):
        for hit in hits_for_read:
            output_bam.write(hit)


def _write_remaining_hits(hits_for_reads, output_bam):
    try:
        while True:
            hits_for_read = hits_for_reads.next()
            _write_hits(hits_for_read, output_bam)
    except StopIteration:
        pass


def _compare_and_write_hits(s1_hits_for_read, s2_hits_for_read,
                            s1_output_bam, s2_output_bam):
    # Compare the hits for a particular read in each species and decide whether
    # the read can be assigned to one species or another, or if it must be
    # rejected as ambiguous
    pass


def _filter_sample_reads(logger, options):
    s1_hits = samutils.open_samfile_for_read(options[SPECIES_ONE_INPUT_BAM])
    s2_hits = samutils.open_samfile_for_read(options[SPECIES_TWO_INPUT_BAM])

    s1_hits_for_reads = _hits_generator(s1_hits)
    s2_hits_for_reads = _hits_generator(s2_hits)

    s1_output_bam = samutils.open_samfile_for_write(
        options[SPECIES_ONE_OUTPUT_BAM], s1_hits)
    s2_output_bam = samutils.open_samfile_for_write(
        options[SPECIES_ONE_OUTPUT_BAM], s2_hits)

    s1_hits_for_read = None
    s2_hits_for_read = None

    while True:
    # 1. Attempt to read all hits for the first/next read in each species'
    # input BAM file. Go to (2).
        if s1_hits_for_read is None:
            try:
                s1_hits_for_read = s1_hits_for_reads.next()
            except StopIteration:
    # 2. If no more reads can be read for species 1:
    #       - all remaining reads in the input file for species 2 can be written
    #       to the output file for species 2, or discarded as inherently
    #       potentially ambiguous.
    #       - END.
    #    Else go to (3).
                _write_remaining_hits(s2_hits_for_read, s2_output_bam)
                break

        if s2_hits_for_read is None:
            try:
                s2_hits_for_read = s2_hits_for_reads.next()
            except StopIteration:
    # 3. If no more reads can be read for species 2:
    #       - all remaining reads in the input file for species 1 can be
    #       written to the output file for species 1, or discarded as
    #       inherently potentially ambiguous.
    #       - END.
    #   Else go to (4).
                _write_remaining_hits(s1_hits_for_read, s1_output_bam)
                break

    # 4. We have the hits for a read for each species. Examine the names of the
    # reads in species 1 and species 2. Go to (5).
        s1_read_name = s1_hits_for_read[0].query_name
        s2_read_name = s2_hits_for_read[0].query_name

    # 5. If the read names are the same:
        if s1_read_name == s2_read_name:
    #       - compare the hits for each species to determine which species to
    #       assign the read to.
    #       - write the hits to the output file for the correct species (or
    #       neither, if the read cannot be unambiguously assigned).
    #       - Go to (1).
    #   Else go to (6).
            _compare_and_write_hits(s1_hits_for_read, s2_hits_for_read,
                                    s1_output_bam, s2_output_bam)
            s1_hits_for_read = None
            s2_hits_for_read = None
    # 6. If species 1 read name < species 2 read name, there are no hits for
    # the species 1 read in species 2:
        elif s1_read_name < s2_read_name:
    #       - write the hits for this read to the output file for species 1, or
    #       discard as inherently potentially ambiguous.
    #       - attempt to read hits for the next read for species 1. If no more
    #       reads can be read for species 1, go to (2). Else go to (4).
    #   Else go to (7)
            _write_hits(s1_hits_for_read, s1_output_bam)
            s1_hits_for_read = None
    # 7. If species 2 read name < species 1 read name, there are no hits for
    # the species 2 read in species 1:
        else:
    #       - write the hits for this read to the output file for species 2, or
    #       discard as inherently potentially ambiguous.
    #       - attempt to read hits for the next read for species 2. If no more
    #       reads can be read for species 2, go to (3). Else go to (4).
            _write_hits(s2_hits_for_read, s2_output_bam)
            s2_hits_for_read = None
    # n.b. throughout, we should always keep a record of the last read name
    # read for a species; when hits for the next read are read for that
    # species, we should check the reads are in name order, and exit with a
    # failure code if not.


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

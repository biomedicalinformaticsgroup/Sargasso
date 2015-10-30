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
import sys

from . import options as opt
from . import samutils
from .__init__ import __version__

SPECIES_ONE = "<species-one>"
SPECIES_ONE_INPUT_BAM = "<species-one-input-bam>"
SPECIES_ONE_OUTPUT_BAM = "<species-one-output-bam>"
SPECIES_TWO = "<species-two>"
SPECIES_TWO_INPUT_BAM = "<species-two-input-bam>"
SPECIES_TWO_OUTPUT_BAM = "<species-two-output-bam>"


ASSIGNED_TO_SPECIES_ONE = 1
ASSIGNED_TO_SPECIES_TWO = 2
ASSIGNED_TO_NEITHER_REJECTED = 3
ASSIGNED_TO_NEITHER_AMBIGUOUS = 4


class _SeparationStats:
    def __init__(self, name):
        self.name = name
        self.hits_written = 0
        self.reads_written = 0
        self.hits_rejected = 0
        self.reads_rejected = 0
        self.hits_ambiguous = 0
        self.reads_ambiguous = 0

    def filtered_hits(self, hits):
        self.hits_written += len(hits)
        self.reads_written += 1

    def rejected_hits(self, hits):
        self.hits_rejected += len(hits)
        self.reads_rejected += 1

    def ambiguous_hits(self, hits):
        self.hits_ambiguous += len(hits)
        self.reads_ambiguous += 1

    def __str__(self):
        return ("{n}: wrote {f} filtered hits for {fr} reads; {r} hits for " +
                "{rr} reads were rejected outright, and {a} hits for " +
                "{ar} reads were rejected as ambiguous.").format(
            n=self.name,
            f=self.hits_written, fr=self.reads_written,
            r=self.hits_rejected, rr=self.reads_rejected,
            a=self.hits_ambiguous, ar=self.reads_ambiguous)


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


class _HitsInfo:
    def __init__(self, hits):
        self.hits = hits
        self.multimaps = None
        self.max_mismatches = None
        self.min_mismatches = None
        self.cigar = None

    def get_multimaps(self):
        if self.multimaps is None:
            self.multimaps = samutils.get_multimaps(self.hits[0])
        return self.multimaps

    def get_max_mismatches(self):
        if self.max_mismatches is None:
            self.max_mismatches = 0
            for hit in self.hits:
                mismatches = samutils.get_mismatches(hit)
                if mismatches > self.max_mismatches:
                    self.max_mismatches = mismatches
        return self.max_mismatches

    def get_min_mismatches(self):
        if self.min_mismatches is None:
            self.min_mismatches = sys.maxsize
            for hit in self.hits:
                mismatches = samutils.get_mismatches(hit)
                if mismatches < self.min_mismatches:
                    self.min_mismatches = mismatches
        return self.min_mismatches

    def get_cigars(self):
        if self.cigar is None:
            self.cigar = {}
            for hit in self.hits:
                self.cigar[hit] = samutils.get_cigar(hit)
        return self.cigar.values()


class _HitsChecker:
    def check_hits(self, hits_info):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species. For example, in the ultra-conservative
        # strategy, this means:
        # - no multi-mapping
        # - no mismatches
        # - CIGAR string satisfactory
        if hits_info.get_multimaps() > 1:
            return False

        if hits_info.get_max_mismatches() > 0:
            return False

        if not self.check_cigar(hits_info):
            return False

        return True

    def assign_hits(self, s1_hits_info, s2_hits_info):
        # TODO: currently making the assumption that min = max mismatches
        # TODO: initially try to get this to behave in the same way as the original code
        s1_multimaps = s1_hits_info.get_multimaps()
        s2_multimaps = s2_hits_info.get_multimaps()

        if s1_multimaps > 1 or s2_multimaps > 1:
            return ASSIGNED_TO_NEITHER_REJECTED

        s1_mismatches = s1_hits_info.get_max_mismatches()
        s2_mismatches = s2_hits_info.get_max_mismatches()

        if s1_mismatches < s2_mismatches:
            return ASSIGNED_TO_SPECIES_ONE
        elif s2_mismatches < s1_mismatches:
            return ASSIGNED_TO_SPECIES_TWO
        else:
            s1_cigar_check = self.check_cigar(s1_hits_info)
            s2_cigar_check = self.check_cigar(s2_hits_info)

            if s1_cigar_check and not s2_cigar_check:
                return ASSIGNED_TO_SPECIES_ONE
            elif not s1_cigar_check and s2_cigar_check:
                return ASSIGNED_TO_SPECIES_TWO

        return ASSIGNED_TO_NEITHER_AMBIGUOUS

    def check_cigar(self, hits_info):
        for cigar in hits_info.get_cigars():
            for c in ["I", "D", "S", "H", "P", "X"]:
                if c in cigar:
                    return False
            if "N" in cigar:
                segments = cigar.split("M")
                if int(segments[0]) < 5 or int(segments[1].split("N")[1]) < 5:
                    return False
        return True


def _hits_generator(all_hits):
    last_hit_name = None
    current_hits = []

    for hit in samutils.all_hits(all_hits):
        if last_hit_name is None:
            last_hit_name = hit.query_name

        if hit.query_name < last_hit_name:
            # TODO: throw an exception if hits are out of read name order, and
            # handle this gracefully
            pass
        elif hit.query_name == last_hit_name:
            current_hits.append(hit)
        else:
            yield current_hits
            last_hit_name = hit.query_name
            current_hits = [hit]

    yield current_hits


def _write_hits(hits_for_read, output_bam):
    for hit in hits_for_read:
        output_bam.write(hit)
        pass


def _write_remaining_hits(hits_checker, hits_for_reads, stats, output_bam):
    try:
        while True:
            hits_for_read = hits_for_reads.next()
            if hits_checker.check_hits(_HitsInfo(hits_for_read)):
                stats.filtered_hits(hits_for_read)
                _write_hits(hits_for_read, output_bam)
            else:
                stats.rejected_hits(hits_for_read)
    except StopIteration:
        pass


def _compare_and_write_hits(s1_hits_for_read, s2_hits_for_read,
                            s1_output_bam, s2_output_bam,
                            s1_stats, s2_stats, hits_checker):
    # Compare the hits for a particular read in each species and decide whether
    # the read can be assigned to one species or another, or if it must be
    # rejected as ambiguous
    s1_hits_info = _HitsInfo(s1_hits_for_read)
    s2_hits_info = _HitsInfo(s2_hits_for_read)

    assignment = hits_checker.assign_hits(s1_hits_info, s2_hits_info)

    if assignment == ASSIGNED_TO_SPECIES_ONE:
        s2_stats.rejected_hits(s2_hits_for_read)
        if hits_checker.check_hits(s1_hits_info):
            s1_stats.filtered_hits(s1_hits_for_read)
            _write_hits(s1_hits_for_read, s1_output_bam)
        else:
            s1_stats.rejected_hits(s1_hits_for_read)
    elif assignment == ASSIGNED_TO_SPECIES_TWO:
        s1_stats.rejected_hits(s1_hits_for_read)
        if hits_checker.check_hits(s2_hits_info):
            s2_stats.filtered_hits(s2_hits_for_read)
            _write_hits(s2_hits_for_read, s2_output_bam)
        else:
            s2_stats.rejected_hits(s2_hits_for_read)
    elif assignment == ASSIGNED_TO_NEITHER_REJECTED:
        s1_stats.rejected_hits(s1_hits_for_read)
        s2_stats.rejected_hits(s2_hits_for_read)
    else:
        s1_stats.ambiguous_hits(s1_hits_for_read)
        s2_stats.ambiguous_hits(s2_hits_for_read)


def _filter_sample_reads(logger, options):
    logger.info("Starting species separation.")

    s1_stats = _SeparationStats("Species 1")
    s2_stats = _SeparationStats("Species 2")

    s1_hits = samutils.open_samfile_for_read(options[SPECIES_ONE_INPUT_BAM])
    s2_hits = samutils.open_samfile_for_read(options[SPECIES_TWO_INPUT_BAM])

    s1_hits_for_reads = _hits_generator(s1_hits)
    s2_hits_for_reads = _hits_generator(s2_hits)

    s1_output_bam = samutils.open_samfile_for_write(
        options[SPECIES_ONE_OUTPUT_BAM], s1_hits)
    s2_output_bam = samutils.open_samfile_for_write(
        options[SPECIES_TWO_OUTPUT_BAM], s2_hits)

    s1_hits_for_read = None
    s2_hits_for_read = None

    s1_count = 0
    s2_count = 0

    hits_checker = _HitsChecker()

    while True:
    # 1. Attempt to read all hits for the first/next read in each species'
    # input BAM file. Go to (2).
        if s1_hits_for_read is None:
            try:
                s1_hits_for_read = s1_hits_for_reads.next()
                s1_count += 1
                if s1_count % 1000000 == 0:
                    logger.debug("Read " + str(s1_count) + " reads from species 1")
            except StopIteration:
    # 2. If no more reads can be read for species 1:
    #       - all remaining reads in the input file for species 2 can be written
    #       to the output file for species 2, or discarded as inherently
    #       potentially ambiguous.
    #       - END.
    #    Else go to (3).
                _write_remaining_hits(hits_checker, s2_hits_for_reads, s2_stats, s2_output_bam)
                break

        if s2_hits_for_read is None:
            try:
                s2_hits_for_read = s2_hits_for_reads.next()
                s2_count += 1
                if s2_count % 1000000 == 0:
                    logger.debug("Read " + str(s2_count) + " reads from species 2")
            except StopIteration:
    # 3. If no more reads can be read for species 2:
    #       - all remaining reads in the input file for species 1 can be
    #       written to the output file for species 1, or discarded as
    #       inherently potentially ambiguous.
    #       - END.
    #   Else go to (4).
                _write_remaining_hits(hits_checker, s1_hits_for_reads, s1_stats, s1_output_bam)
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
                                    s1_output_bam, s2_output_bam,
                                    s1_stats, s2_stats, hits_checker)
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
            if hits_checker.check_hits(_HitsInfo(s1_hits_for_read)):
                s1_stats.filtered_hits(s1_hits_for_read)
                _write_hits(s1_hits_for_read, s1_output_bam)
            else:
                s1_stats.rejected_hits(s1_hits_for_read)
            s1_hits_for_read = None
    # 7. If species 2 read name < species 1 read name, there are no hits for
    # the species 2 read in species 1:
        else:
    #       - write the hits for this read to the output file for species 2, or
    #       discard as inherently potentially ambiguous.
    #       - attempt to read hits for the next read for species 2. If no more
    #       reads can be read for species 2, go to (3). Else go to (4).
            if hits_checker.check_hits(_HitsInfo(s2_hits_for_read)):
                s2_stats.filtered_hits(s2_hits_for_read)
                _write_hits(s2_hits_for_read, s2_output_bam)
            else:
                s2_stats.rejected_hits(s2_hits_for_read)
            s2_hits_for_read = None
    # n.b. throughout, we should always keep a record of the last read name
    # read for a species; when hits for the next read are read for that
    # species, we should check the reads are in name order, and exit with a
    # failure code if not.

    logger.info(s1_stats)
    logger.info(s2_stats)


def filter_reads_parallel(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="filter_sample_reads v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    _filter_sample_reads(logger, options)

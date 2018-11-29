import os
import os.path
import os.path
import schema
import sargasso.separator.options as opts

from sargasso.filter import filterer, hits_checker
from sargasso.separator.commandline_parser import CommandlineParser
from sargasso.separator.parameter_validator import ParameterValidator
from sargasso.utils import log


class SampleFilter(object):
    DOC = """Usage:
    filter_sample_reads -h | --help
    filter_sample_reads -v | --version
    filter_sample_reads <data-type> [<args>...]

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}

The available commands are:
   rnaseq       RNA-seq data
   chipseq      ChIP-seq data

"""
    SPECIES_INPUT_BAM = "<species-input-bam>"
    SPECIES_OUTPUT_BAM = "<species-output-bam>"

    def __init__(self, filterer_cls, hits_checker_cls, commandline_parser):

        self.filterer_cls = filterer_cls
        self.hits_checker_cls = hits_checker_cls
        self.commandline_parser = commandline_parser

    @classmethod
    def _validate_command_line_options(cls, options):
        try:
            ParameterValidator.validate_log_level(options)
            for index, species in enumerate(options[opts.SPECIES_ARG]):
                ParameterValidator.validate_file_option(
                    options[SampleFilter.SPECIES_INPUT_BAM][index],
                    "Could not find input BAM file for species {i}".format(i=index))

            ParameterValidator.validate_threshold_options(
                options,
                opts.MISMATCH_THRESHOLD_ARG,
                opts.MINMATCH_THRESHOLD_ARG,
                opts.MULTIMAP_THRESHOLD_ARG)

        except schema.SchemaError as exc:
            exit(exc.code)

    def _filter_sample_reads(self, logger, options):
        logger.info("Starting species separation.")

        h_check = self.hits_checker_cls(
                options[opts.MISMATCH_THRESHOLD_ARG],
                options[opts.MINMATCH_THRESHOLD_ARG],
                options[opts.MULTIMAP_THRESHOLD_ARG],
                options[opts.REJECT_MULTIMAPS],
                logger)

        filterers = [self.filterer_cls(i + 1,
                                       options[SampleFilter.SPECIES_INPUT_BAM][i],
                                       options[SampleFilter.SPECIES_OUTPUT_BAM][i],
                                       logger)
                     for i, species in enumerate(options[opts.SPECIES_ARG])]

        all_filterers = filterers

        while True:
            # Retain only filterers which have hits for the remaining reads
            filterers = [f for f in filterers if self._get_next_read_name(f) is not None]

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
            min_read_name = self._get_next_read_name(filterers[0])

            for cfilt in filterers[1:]:
                read_name = self._get_next_read_name(cfilt)
                if read_name == min_read_name:
                    competing_filterers.append(cfilt)
                elif read_name < min_read_name:
                    competing_filterers = [cfilt]
                    min_read_name = read_name

            if __debug__:
                logger.debug("Read:{}".format(competing_filterers[0].hits_for_read[0].qname))

            # If there's only one filterer for this read, write hits for that read
            # to the output file for that species (or discard as ambiguous)
            if len(competing_filterers) == 1:
                if __debug__:
                    logger.debug('assigned due to only one competing filterer!')
                h_check.check_and_write_hits_for_read(competing_filterers[0])
                continue

            # Otherwise compare the hits for each species to determine which
            # species to assign the read to.
            h_check.compare_and_write_hits(competing_filterers)

        for filt in all_filterers:
            filt.log_stats()

        self._write_stats(all_filterers, options[SampleFilter.SPECIES_OUTPUT_BAM][0])

    # write filter stats to table in file
    @classmethod
    def _write_stats(cls, filterers, out_bam):

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

    @classmethod
    def _get_next_read_name(cls, f):
        read_name = None

        try:
            read_name = f.get_next_read_name()
        except StopIteration:
            pass

        return read_name

    def run(self, args):
        # Read in command-line options
        options = self.commandline_parser.parse(args, self.DOC)

        # Validate command-line options
        self._validate_command_line_options(options)

        # Set up logger
        self.logger = log.get_logger_for_options(options)

        self._filter_sample_reads(self.logger, options)


class RnaseqSampleFilter(SampleFilter):
    DOC = """
Usage:
filter_sample_reads <data-type>
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

    def __init__(self, commandline_parser):
        SampleFilter.__init__(
            self, filterer.RnaseqFilterer, hits_checker.RnaseqHitsChecker,
            commandline_parser)


class ChipseqSampleFilter(SampleFilter):
    DOC = """
Usage:
filter_sample_reads <data-type>
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

    def __init__(self, commandline_parser):
        SampleFilter.__init__(
            self, filterer.ChipseqFilterer, hits_checker.ChipseqHitChecker,
            commandline_parser)

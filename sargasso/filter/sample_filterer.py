import os
import os.path
import os.path
import schema

import sargasso.separator.options as opts
from sargasso.filter import hits_manager, hits_checker
from sargasso.separator.parameter_validator import ParameterValidator
from sargasso.utils import log


class SampleFilterer(object):
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
   rnaseq       RNA-sequencing data
   dnaseq       DNA-sequencing data

"""
    SPECIES_INPUT_BAM = "<species-input-bam>"
    SPECIES_OUTPUT_BAM = "<species-output-bam>"

    def __init__(self, hits_manager_cls, hits_checker_cls, commandline_parser):
        self.hits_manager_cls = hits_manager_cls
        self.hits_checker_cls = hits_checker_cls
        self.commandline_parser = commandline_parser

    @classmethod
    def _validate_command_line_options(cls, options):
        try:
            ParameterValidator.validate_log_level(options)
            for index, species in enumerate(options[opts.SPECIES_ARG]):
                ParameterValidator.validate_file_option(
                    options[SampleFilterer.SPECIES_INPUT_BAM][index],
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

        if __debug__:
            outlog = os.path.join( os.path.dirname(options[SampleFilterer.SPECIES_OUTPUT_BAM][0]),
                                   os.path.basename(options[SampleFilterer.SPECIES_OUTPUT_BAM][0]).split('___')[-2] + '.log')
            fhandler = log.get_debug_file_handler(outlog)
            logger.addHandler(fhandler)

        h_check = self.hits_checker_cls(
            options[opts.MISMATCH_THRESHOLD_ARG],
            options[opts.MINMATCH_THRESHOLD_ARG],
            options[opts.MULTIMAP_THRESHOLD_ARG],
            options[opts.REJECT_MULTIMAPS],
            logger)

        hits_managers = [self.hits_manager_cls(
            i + 1,
            options[SampleFilterer.SPECIES_INPUT_BAM][i],
            options[SampleFilterer.SPECIES_OUTPUT_BAM][i],
            logger)
            for i, species in enumerate(options[opts.SPECIES_ARG])]

        all_hits_managers = hits_managers

        while True:
            # Retain only hits managers which have hits for the remaining reads
            hits_managers = [m for m in hits_managers
                             if self._get_next_read_name(m) is not None]

            # If no hits managers remain, we're done
            if len(hits_managers) == 0:
                break


            # If only one hits manager remains, all remaining reads in the
            # input file for that species can be written to the output file for
            # that species (or discarded as ambiguous, if necessary).
            if len(hits_managers) == 1:
                ## xintodo debug remove
                ## to ignore one hits_manager case
                # hits_managers[0].clear_hits()
                # continue
                ## only look into specific read
                # if hits_managers[0].hits_for_read[0].qname != 'SRR5467502.10000_10000_length=100':
                #    hits_managers[0].clear_hits()
                # continue
                h_check.check_and_write_hits_for_remaining_reads(hits_managers[0])
                break

            # Compare read names for the hits managers, and find the set that
            # have the "lowest" name
            competing_hits_managers = [hits_managers[0]]
            min_read_name = self._get_next_read_name(hits_managers[0])

            for cman in hits_managers[1:]:
                read_name = self._get_next_read_name(cman)
                if read_name == min_read_name:
                    competing_hits_managers.append(cman)
                elif read_name < min_read_name:
                    competing_hits_managers = [cman]
                    min_read_name = read_name

            ## xintodo debug remove
            ## only look into specific read
            # if min_read_name != 'SRR5467502.28463_28463_length=100':
            #     for hits_manager in competing_hits_managers:
            #         hits_manager.clear_hits()
            #     continue
            # only look into multiple competing_hits_managers
            # if len(competing_hits_managers) < 2:
            #     for hits_manager in competing_hits_managers:
            #         hits_manager.clear_hits()
            #     continue

            if __debug__: logger.debug("({}) Read:{}, Species:{}".format(
                len(competing_hits_managers),
                competing_hits_managers[0].hits_for_read[0].qname,
                [m.species_id for m in competing_hits_managers]))


            # If there's only one hits manager for this read, write hits for
            # that read to the output file for that species (or discard as
            # ambiguous)
            if len(competing_hits_managers) == 1:
                h_check.check_and_write_hits_for_read(competing_hits_managers[0])
                continue

            # Otherwise compare the hits for each species to determine which
            # species to assign the read to.
            h_check.compare_and_write_hits(competing_hits_managers)

        for filt in all_hits_managers:
            filt.log_stats()

        self._write_stats(all_hits_managers,
                          options[SampleFilterer.SPECIES_OUTPUT_BAM][0])

    # write filter stats to table in file
    @classmethod
    def _write_stats(cls, hits_managers, out_bam):

        stats = []

        for man in hits_managers:
            mstats = man.stats

            stats += [mstats.hits_written, mstats.reads_written,
                      mstats.hits_rejected, mstats.reads_rejected,
                      mstats.hits_ambiguous, mstats.reads_ambiguous]

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


class RnaSeqSampleFilterer(SampleFilterer):
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
of mixed species sequencing reads against a number of species' genomes, and
determines, if possible, from which species each read originates. Disambiguated
read mappings are written to species-specific output BAM files.

In normal operation, the user should not need to execute this script by hand
themselves.

Note: the input BAM files MUST be sorted in read name order. Failure to
ensure input BAM files are correctly sorted will result in erroneous output.
"""

    def __init__(self, commandline_parser):
        SampleFilterer.__init__(
            self, hits_manager.RnaSeqHitsManager, hits_checker.RnaSeqHitsChecker, commandline_parser)


class DnaSeqSampleFilterer(SampleFilterer):
    DOC = RnaSeqSampleFilterer.DOC

    def __init__(self, commandline_parser):
        SampleFilterer.__init__(
            self, hits_manager.DnaSeqHitsManager, hits_checker.DnaSeqHitsChecker, commandline_parser)


class BisulfiteSampleFilterer(SampleFilterer):
    DOC = RnaSeqSampleFilterer.DOC

    def __init__(self, commandline_parser):
        SampleFilterer.__init__(
            self, hits_manager.BisulfiteHitsManager, hits_checker.BisulfiteHitsChecker, commandline_parser)

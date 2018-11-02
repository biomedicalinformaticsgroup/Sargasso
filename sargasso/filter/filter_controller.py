import subprocess

import os
import os.path
import schema

from sargasso.separator.commandline_parser import CommandlineParser
from sargasso.separator.commandline_parser import CommandlineParserManager
from sargasso.separator.options import Options
from sargasso.separator.parameter_validator import ParameterValidator
from sargasso.utlis.factory import Manager
from sargasso.utlis.log import LoggerManager


class FilterController(object):
    DOC = """
Usage: filter_control [-h] <data-type> [<args>...]
    
options:
   -h
   
The available commands are:
   rnaseq   rnaseq data
   chipseq  chipseq data
"""
    data_type = None

    DATA_TYPE = "<data-type>"
    BLOCK_DIR = "<block-dir>"
    OUTPUT_DIR = "<output-dir>"
    SAMPLE_NAME = "<sample-name>"
    SPECIES = "<species>"
    MISMATCH_THRESHOLD = "<mismatch_threshold>"
    MINMATCH_THRESHOLD = "<minmatch_threshold>"
    MULTIMAP_THRESHOLD = "<multimap_threshold>"
    REJECT_MULTIMAPS = "--reject-multimaps"
    BLOCK_FILE_SEPARATOR = "___"

    def __init__(self, data_type, commandline_parser, logger):
        self.data_type = data_type
        self.commandline_parser = commandline_parser
        self.logger = logger

    def _validate_command_line_options(self, options):
        """
        Validate command line options are correctly specified.

        options: dictionary of command line options.
        """
        try:
            ParameterValidator.validate_log_level(options)
            ParameterValidator.validate_dir_option(
                options[FilterController.BLOCK_DIR],
                "Mapped reads block file directory does not exist")
            ParameterValidator.validate_dir_option(
                options[FilterController.OUTPUT_DIR],
                "Filtered reads output directory does not exist")
        except schema.SchemaError as exc:
            exit("Exiting: " + exc.code)

    def _all_processes_finished(self, processes):
        """
        Return True if all processes have finished.

        processes: List of subprocess instances.
        """
        finished = True
        i = 0

        while True:
            if i < len(processes):
                if processes[i].poll() is None:
                    finished = False
                elif processes[i].poll() == 0:
                    del processes[i]
                    i = i - 1
                if not finished:
                    return False
                i = i + 1
            else:
                break

        return True

    @staticmethod
    def _get_block_files(block_dir, sp1):
        is_block_file = lambda f: not os.path.isdir(os.path.join(block_dir, f))
        block_files = [f for f in os.listdir(block_dir) if is_block_file(f)]
        block_files_out = []

        for block_file in sorted(block_files):
            sections = block_file.split(FilterController.BLOCK_FILE_SEPARATOR)
            if sections[1] == sp1:
                block_files_out.append(block_file)

        return block_files_out

    @staticmethod
    def _get_input_path(block_dir, block_file, species):
        sections = block_file.split(FilterController.BLOCK_FILE_SEPARATOR)
        sections[1] = species
        block_file = FilterController.BLOCK_FILE_SEPARATOR.join(sections)
        return os.path.join(block_dir, block_file)

    @staticmethod
    def _initialise_result_file(options):
        """
        Initialise results summary file.

        out_dir: Directory into which filtered BAM files will be written.
        """
        cols = []

        for index, species in enumerate(options[FilterController.SPECIES]):
            species_text = str(species)

            cols += [
                "Assigned-Hits-" + species_text, "Assigned-Reads-" + species_text,
                "Rejected-Hits-" + species_text, "Rejected-Reads-" + species_text,
                "Ambiguous-Hits-" + species_text, "Ambiguous-Reads-" + species_text
            ]

        out_file = os.path.join(options[FilterController.OUTPUT_DIR], "filtering_result_summary.txt")
        with open(out_file, 'w') as outf:
            outf.write("\t".join(cols) + "\n")

    def _run_processes(self, logger, options):
        """
        Run filtering script in a separate process for each pair of block files.

        logger: logging object
        options: dictionary of command-line options
        """
        block_files = self._get_block_files(options[FilterController.BLOCK_DIR], options[FilterController.SPECIES][0])

        # keep track of all processes
        all_processes = []
        proc_no = -1

        # initialise results file
        self._initialise_result_file(options)

        # cycle through chunks
        for block_file in block_files:
            proc_no += 1

            commands = ["filter_sample_reads", self.data_type,
                        options[FilterController.MISMATCH_THRESHOLD], options[FilterController.MINMATCH_THRESHOLD],
                        options[FilterController.MULTIMAP_THRESHOLD]]

            for species in options[FilterController.SPECIES]:
                sp_in = self._get_input_path(options[FilterController.BLOCK_DIR], block_file, species)

                get_output_path = lambda x: os.path.join(
                    options[FilterController.OUTPUT_DIR],
                    "___".join([options[FilterController.SAMPLE_NAME], x,
                              str(proc_no), "filtered.bam"]))

                sp_out = get_output_path(species)

                commands += [species, sp_in, os.path.abspath(sp_out)]

            if options[FilterController.REJECT_MULTIMAPS]:
                commands.append("--reject-multimaps")

            proc = subprocess.Popen(commands)
            all_processes.append(proc)

        # check all processes finished
        while not self._all_processes_finished(all_processes):
            logger.info("Waiting for Threads")
            all_processes[0].wait()

        logger.info("Filtering Complete")

    def run(self, args):
        """
        Reads mapped read block files and parallelises their filtering.

        args: list of command line arguments
        """
        # Read in command line options
        options = self.commandline_parser.parse(args, self.DOC)

        # Validate command line options
        self._validate_command_line_options(options)

        # Set up logger
        self.logger.init(options[Options.LOG_LEVEL_OPTION])

        # Parallelise species separation filtering of mapped read block files
        self._run_processes(self.logger, options)


class RnaseqFilterController(FilterController):
    DOC = """
Usage:
    filter_control <data-type>
        [--log-level=<log-level>] [--reject-multimaps]
        <block-dir> <output-dir> <sample-name>
        <mismatch_threshold> <minmatch_threshold> <multimap_threshold>
        (<species>) (<species>) ...

Option:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<block-dir>
    Directory containing pairs of mapped read BAM files.
<output-dir>
    Directory into which species-separated reads will be written.
<sample-name>
    Name of RNA-seq sample being processed.
<species>
    Name of species.
<mismatch-threshold>
    Maximum percentage of read bases allowed to be mismatches against the
    genome during filtering.
<minmatch-threshold>
    Maximum percentage of read length allowed to not be mapped during
    filtering.
<multimap_threshold>
    Maximum number of multiple mappings allowed during filtering.
--reject-multimaps
    If set, any read which multimaps to any species' genome will be rejected
    and not be assigned to any species.

filter_control takes a directory containing sets of BAM files as input, each
set being the result of mapping a set of mixed species RNA-seq reads against
a number of species' genomes. Each set of BAM files is passed to an instance of
the script filter_sample_reads, running on a separate thread, which determines
where possible from which species each read originates. Read mappings for each
pair of input files are written to a set of species-specific output BAM files
in the specified output directory.

In normal operation, the user should not need to execute this script by hand
themselves.

Note: the input BAM files MUST be sorted in read name order. Failure to ensure
input BAM files are correctly sorted will result in erroneous output.
"""
    pass


class ChipseqFilterController(FilterController):
    DOC = """
Usage:
    filter_control <data-type>
        [--log-level=<log-level>] [--reject-multimaps]
        <block-dir> <output-dir> <sample-name>
        <mismatch_threshold> <minmatch_threshold> <multimap_threshold>
        (<species>) (<species>) ...

Option:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<block-dir>
    Directory containing pairs of mapped read BAM files.
<output-dir>
    Directory into which species-separated reads will be written.
<sample-name>
    Name of RNA-seq sample being processed.
<species>
    Name of species.
<mismatch-threshold>
    Maximum percentage of read bases allowed to be mismatches against the
    genome during filtering.
<minmatch-threshold>
    Maximum percentage of read length allowed to not be mapped during
    filtering.
<multimap_threshold>
    Maximum number of multiple mappings allowed during filtering.
--reject-multimaps
    If set, any read which multimaps to any species' genome will be rejected
    and not be assigned to any species.

filter_control takes a directory containing sets of BAM files as input, each
set being the result of mapping a set of mixed species RNA-seq reads against
a number of species' genomes. Each set of BAM files is passed to an instance of
the script filter_sample_reads, running on a separate thread, which determines
where possible from which species each read originates. Read mappings for each
pair of input files are written to a set of species-specific output BAM files
in the specified output directory.

In normal operation, the user should not need to execute this script by hand
themselves.

Note: the input BAM files MUST be sorted in read name order. Failure to ensure
input BAM files are correctly sorted will result in erroneous output.
"""
    pass


class FilterControllerManager(Manager):
    FILTERCONTROLLER = {"rnaseq": RnaseqFilterController,
                        "chipseq": ChipseqFilterController}

    def __init__(self):
        pass

    @classmethod
    def _create(cls, data_type):
        commandline_parser = CommandlineParserManager.get(data_type)
        logger = LoggerManager.get()
        return cls.FILTERCONTROLLER[data_type](data_type, commandline_parser, logger)

    @classmethod
    def get(cls, data_type):
        return cls._create(data_type)


def filter_control(args):
    data_type = CommandlineParser.parse_datatype(args, FilterController.DOC, FilterController.DATA_TYPE)
    try:
        ParameterValidator.validate_datatype(data_type)
    except schema.SchemaError as exc:
        exit("Exiting: " + exc.code)

    filter_controller = FilterControllerManager.get(data_type)
    filter_controller.run(args)

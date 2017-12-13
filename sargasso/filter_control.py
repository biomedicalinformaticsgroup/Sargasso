#!/usr/bin/env python

"""Usage:
    filter_control
        [--log-level=<log-level>] [--reject-multimaps]
        <block-dir> <output-dir> <sample-name>
        <mismatch_threshold> <minmatch_threshold>
        <multimap_threshold> <overhang-threshold>
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
<overhang-threshold>
    The minimum number of bases that are allowed on
    either side of an exon boundary for a read mapping to be accepted
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

import docopt
import os
import os.path
import schema
import subprocess

from . import options as opt
from .__init__ import __version__

BLOCK_DIR = "<block-dir>"
OUTPUT_DIR = "<output-dir>"
SAMPLE_NAME = "<sample-name>"
SPECIES = "<species>"
MISMATCH_THRESHOLD = "<mismatch_threshold>"
MINMATCH_THRESHOLD = "<minmatch_threshold>"
MULTIMAP_THRESHOLD = "<multimap_threshold>"
REJECT_MULTIMAPS = "--reject-multimaps"
OVERHANG_THRESHOLD = "<overhang-threshold>"

BLOCK_FILE_SEPARATOR = "___"


def _validate_command_line_options(options):
    """
    Validate command line options are correctly specified.

    options: dictionary of command line options.
    """
    try:
        opt.validate_log_level(options)
        opt.validate_dir_option(
            options[BLOCK_DIR],
            "Mapped reads block file directory does not exist")
        opt.validate_dir_option(
            options[OUTPUT_DIR],
            "Filtered reads output directory does not exist")
    except schema.SchemaError as exc:
        exit("Exiting: " + exc.code)


def _all_processes_finished(processes):
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


def _get_block_files(block_dir, sp1):
    is_block_file = lambda f: not os.path.isdir(os.path.join(block_dir, f))
    block_files = [f for f in os.listdir(block_dir) if is_block_file(f)]
    block_files_out = []

    for block_file in sorted(block_files):
        sections = block_file.split(BLOCK_FILE_SEPARATOR)
        if sections[1] == sp1:
            block_files_out.append(block_file)

    return block_files_out


def _get_input_path(block_dir, block_file, species):
    sections = block_file.split(BLOCK_FILE_SEPARATOR)
    sections[1] = species
    block_file = BLOCK_FILE_SEPARATOR.join(sections)

    return os.path.join(block_dir, block_file)


def _run_processes(logger, options):
    """
    Run filtering script in a separate process for each pair of block files.

    logger: logging object
    options: dictionary of command-line options
    """
    block_files = _get_block_files(options[BLOCK_DIR], options[SPECIES][0])

    # keep track of all processes
    all_processes = []
    proc_no = -1

    # initialise results file
    _initialise_result_file(options)

    # cycle through chunks
    for block_file in block_files:
        proc_no += 1

        commands = ["filter_sample_reads",
                    options[MISMATCH_THRESHOLD], options[MINMATCH_THRESHOLD],
                    options[MULTIMAP_THRESHOLD], options[OVERHANG_THRESHOLD]]

        for species in options[SPECIES]:
            sp_in = _get_input_path(options[BLOCK_DIR], block_file, species)

            get_output_path = lambda x: os.path.join(
                options[OUTPUT_DIR],
                "_".join([options[SAMPLE_NAME], x,
                         str(proc_no), "filtered.bam"]))

            sp_out = get_output_path(species)

            commands += [species, sp_in, os.path.abspath(sp_out)]

        if options[REJECT_MULTIMAPS]:
            commands.append("--reject-multimaps")

        proc = subprocess.Popen(commands)
        all_processes.append(proc)

    # check all processes finished
    while not _all_processes_finished(all_processes):
        logger.info("Waiting for Threads")
        all_processes[0].wait()

    logger.info("Filtering Complete")


def _initialise_result_file(options):
    """
    Initialise results summary file.

    out_dir: Directory into which filtered BAM files will be written.
    """
    cols = []

    for index, species in enumerate(options[SPECIES]):
        species_text = str(species)

        cols += [
            "Assigned-Hits-" + species_text, "Assigned-Reads-" + species_text,
            "Rejected-Hits-" + species_text, "Rejected-Reads-" + species_text,
            "Ambiguous-Hits-" + species_text, "Ambiguous-Reads-" + species_text
        ]

    out_file = os.path.join(options[OUTPUT_DIR], "filtering_result_summary.txt")
    with open(out_file, 'w') as outf:
        outf.write("\t".join(cols) + "\n")


def filter_control(args):
    """
    Reads mapped read block files and parallelises their filtering.

    args: list of command line arguments
    """
    # Read in command line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="filter_control v" + __version__)

    # Validate command line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Parallelise species separation filtering of mapped read block files
    _run_processes(logger, options)

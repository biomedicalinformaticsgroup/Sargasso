#!/usr/bin/env python

"""Usage:
    filter_control
        [--log-level=<log-level>] [--reject-multimaps]
        <block-dir> <output-dir> <sample-name> <species-one> <species-two>
        <mismatch_threshold> <minmatch_threshold> <multimap_threshold> <overhang-threshold>

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
<species-one>
    Name of first species.
<species-two>
    Name of second species.
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
    If set, any read which multimaps to either species' genome will be rejected
    and not be assigned to either species.

filter_control takes a directory containing pairs of BAM files as input, each
pair being the result of mapping a set of mixed species RNA-seq reads against
the two species' genomes. Each pair of BAM files is passed to an instance of
the script filter_sample_reads, running on a separate thread, which determines
where possible from which species each read originates. Read mappings for each
pair of input files are written to a pair of species-specific output BAM files
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
SPECIES_ONE = "<species-one>"
SPECIES_TWO = "<species-two>"
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
        # TODO: format exit message for 80 columns
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


def _get_block_file_dictionary(block_files, sp1, sp2):
    """
    Return a dictionary mapping from species 1 to species 2 block files.

    block_files: List of block file names
    sp1: species one name
    sp2: species two name
    """
    block_file_dict = {}

    for block_file in sorted(block_files):
	#print "BLOCK: "+str(block_file)
	#print "ALL BLOCKS: "+str(block_files)
        sections = block_file.split(BLOCK_FILE_SEPARATOR)
        if sections[1] == sp1:
            sections[1] = sp2
            block_file_dict[block_file] = BLOCK_FILE_SEPARATOR.join(sections)

    return block_file_dict


def _run_processes(logger, options):
    """
    Run filtering script in a separate process for each pair of block files.

    logger: logging object
    options: dictionary of command-line options
    """
    # get block files
    is_block_file = lambda f: not os.path.isdir(
        os.path.join(options[BLOCK_DIR], f))
    block_files = [f for f in os.listdir(options[BLOCK_DIR])
                   if is_block_file(f)]

    # create dictionary mapping from species 1 to species 2 block files
    block_file_dict = _get_block_file_dictionary(
        block_files, options[SPECIES_ONE], options[SPECIES_TWO])

    # keep track of all processes
    all_processes = []
    proc_no = -1

    get_input_path = lambda x: os.path.join(options[BLOCK_DIR], x)

    # initialise results file
    write_result_file(options[OUTPUT_DIR])

    # cycle through chunks
    for sp1_file in block_file_dict.keys():
        proc_no += 1

        # create input paths
        sp1_in = get_input_path(sp1_file)
        sp2_in = get_input_path(block_file_dict[sp1_file])

        # create output paths
        get_output_path = lambda x: os.path.join(
            options[OUTPUT_DIR],
            "_".join([options[SAMPLE_NAME], x, str(proc_no), "filtered.bam"]))
        sp1_out = get_output_path(options[SPECIES_ONE])
        sp2_out = get_output_path(options[SPECIES_TWO])

        commands = ["filter_sample_reads",
                    options[SPECIES_ONE], sp1_in, os.path.abspath(sp1_out),
                    options[SPECIES_TWO], sp2_in, os.path.abspath(sp2_out),
                    options[MISMATCH_THRESHOLD], options[MINMATCH_THRESHOLD],
                    options[MULTIMAP_THRESHOLD], options[OVERHANG_THRESHOLD]]

        if options[REJECT_MULTIMAPS]:
            commands.append("--reject-multimaps")

        proc = subprocess.Popen(commands)
        all_processes.append(proc)

    # check all processes finished
    while not _all_processes_finished(all_processes):
        logger.info("Waiting for Threads")
        all_processes[0].wait()

    # TODO: need to concatenate output files
    logger.info("Filtering Complete")


# initialise the results file so the threads can append
# DEBUGGING PURPOSES ONLY - Remove once filtering optimisation has been
# completed
def write_result_file(out_dir):
    cols = [
        "Filtered-Hits-S1", "Filtered-Reads-S1",
        "Rejected-Hits-S1", "Rejected-Reads-S1",
        "Ambiguous-Hits-S1", "Ambiguous-Reads-S1",
        "Filtered-Hits-S2", "Filtered-Reads-S2",
        "Rejected-Hits-S2", "Rejected-Reads-S2",
        "Ambiguous-Hits-S2", "Ambiguous-Reads-S2"
    ]

    out_file = os.path.join(out_dir, "filtering_result_summary.txt")
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

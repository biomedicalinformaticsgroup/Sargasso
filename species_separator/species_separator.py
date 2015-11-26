#!/usr/bin/env python

"""Usage:
    species_separator [--log-level=<log-level>] [--reads-base-dir=<reads-base-dir>] [--num-threads=<num-threads>] [--s1-gtf=<species-one-gtf-file>] [--s2-gtf=<species-two-gtf-file>] [--s1-genome-fasta=<species-one-genome-fasta>] [--s2-genome-fasta=<species-two-genome-fasta>] [--s1-index=<species-one-star-index>] [--s2-index=<species-two-star-index>] [--run-separation] <species-one> <species-two> <samples-file> <output-dir>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<species-one>                                   Name of first species.
<species-two>                                   Name of second species.
<samples-file>                                  TSV file giving paths of raw RNA-seq read data files for each sample.
<output-dir>                                    Output directory in which species separation will be performed.
--reads-base-dir=<reads-base-dir>               Base directory for raw RNA-seq read data files.
-t <num-threads> --num-threads=<num-threads>    Number of threads to use for parallel processing. [default: 1]
--s1-gtf=<species-one-gtf-file>                 GTF annotation file for first species.
--s2-gtf=<species-two-gtf-file>                 GTF annotation file for second species.
--s1-genome-fasta=<species-one-genome-fasta>    Directory containing genome FASTA files for first species.
--s2-genome-fasta=<species-two-genome-fasta>    Directory containing genome FASTA files for second species.
--s1-index=<species-one-star-index>             STAR index directory for first species.
--s2-index=<species-two-star-index>             STAR index directory for second species.
--run-separation                                If specified, species separation will be run; otherwise scripts to perform separation will be created but not run.

TODO: what does this script do...
"""

import docopt
import os
import schema

from . import file_writer as fw
from . import options as opt
from .__init__ import __version__

SPECIES_ONE = "<species-one>"
SPECIES_TWO = "<species-two>"
SAMPLES_FILE = "<samples-file>"
OUTPUT_DIR = "<output-dir>"
READS_BASE_DIR = "--reads-base-dir"
NUM_THREADS = "--num-threads"
SPECIES_ONE_GTF = "--s1-gtf"
SPECIES_TWO_GTF = "--s2-gtf"
SPECIES_ONE_GENOME_FASTA = "--s1-genome-fasta"
SPECIES_TWO_GENOME_FASTA = "--s2-genome-fasta"
SPECIES_ONE_INDEX = "--s1-index"
SPECIES_TWO_INDEX = "--s2-index"
RUN_SEPARATION = "--run-separation"

GTF_FILE = "gtf-file"
GENOME_FASTA = "genome-fasta"
STAR_INDEX = "star-index"


def _get_species_options(options, gtf_file_option,
                         genome_fasta_option, star_index_option):
    """
    Return a dictionary containing command-line options for a particular
    species.

    options: dictionary containing all command-line options.
    gtf_file_option: name of option specifying GTF annotation file for the
    species.
    genome_fasta_option: name of option specifying directory containing genome
    FASTA files for the species.
    star_index_option: name of option specifying STAR index directory for the
    species.
    """
    species_options = {}
    species_options[GTF_FILE] = options[gtf_file_option]
    species_options[GENOME_FASTA] = options[genome_fasta_option]
    species_options[STAR_INDEX] = options[star_index_option]
    return species_options


def _validate_species_options(species, species_options):
    """
    Validate command-line options for a particular species are correctly
    specified.

    species: species identification string
    species_options: dictionary of options specific to a particular species.
    """
    opt.validate_file_option(
        species_options[GTF_FILE],
        "Could not open species {species} GTF file".format(species=species),
        nullable=True)
    opt.validate_dir_option(
        species_options[GENOME_FASTA],
        "Genome FASTA directory for species {species} should exist".
        format(species=species),
        nullable=True)
    opt.validate_dir_option(
        species_options[STAR_INDEX],
        "STAR index directory for species {species} should exist".
        format(species=species),
        nullable=True)

    if (species_options[GTF_FILE] is None) != \
            (species_options[GENOME_FASTA] is None):
        raise schema.SchemaError(None, "Should specify both GTF file and " +
                                 "genome FASTA directory for species " +
                                 "{species}.".format(species=species))

    if (species_options[GTF_FILE] is None) == \
            (species_options[STAR_INDEX] is None):
        raise schema.SchemaError(None, "Should specify either GTF file " +
                                 "or STAR index directory for species " +
                                 "{species} (but not both).".
                                 format(species=species))


def _validate_command_line_options(options):
    """
    Validate command line options are correctly specified.

    options: dictionary of command-line options.
    """
    try:
        opt.validate_log_level(options)
        opt.validate_dir_option(
            options[READS_BASE_DIR], "Reads base directory does not exist",
            nullable=True)
        options[NUM_THREADS] = opt.validate_int_option(
            options[NUM_THREADS],
            "Number of threads must be a positive integer",
            min_val=1, nullable=True)
        opt.validate_file_option(
            options[SAMPLES_FILE], "Could not open samples definition file")
        opt.validate_dir_option(
            options[OUTPUT_DIR], "Output directory should not exist",
            should_exist=False)

        species_one_options = _get_species_options(
            options, SPECIES_ONE_GTF, SPECIES_ONE_GENOME_FASTA,
            SPECIES_ONE_INDEX)
        species_two_options = _get_species_options(
            options, SPECIES_TWO_GTF, SPECIES_TWO_GENOME_FASTA,
            SPECIES_TWO_INDEX)

        _validate_species_options("one", species_one_options)
        _validate_species_options("two", species_two_options)

    except schema.SchemaError as exc:
        exit(exc.code)


def _write_variable_definitions(logger, writer, options):
    # TODO: Write variable definitions
    writer.add_line("NUM_THREADS={num_threads}".format(
        num_threads=options[NUM_THREADS]))


def _write_target_variable_definitions(logger, writer, options):
    # TODO: Write target variable definitions
    pass


def _write_phony_targets(logger, writer, options):
    # TODO: Write .PHONY targets
    pass


def _write_all_target(logger, writer, options):
    # TODO: Write "all" target
    pass


def _write_filtered_reads_target(logger, writer, options):
    # TODO: Write "filtered_reads" target
    pass


def _write_sorted_reads_target(logger, writer, options):
    # TODO: Write "sorted_reads" target
    pass


def _write_mapped_reads_target(logger, writer, options):
    # TODO: Write "mapped_reads" target
    pass


def _write_masked_reads_target(logger, writer, options):
    # TODO: Write "masked_reads" target
    pass


def _write_collate_raw_reads_target(logger, writer, options):
    # TODO: Write "collate_raw_reads" target
    pass


def _write_mask_star_index_targets(logger, writer, options):
    # TODO: Write targets for STAR indices for mask sequences
    pass


def _write_main_star_index_targets(logger, writer, options):
    # TODO: Write targets to link to main STAR indices (if we're using
    # pre-built indices) or build main STAR indices (if not)
    pass


def _write_clean_target(logger, writer, options):
    # TODO: Write "clean" target
    pass


def _write_makefile(logger, options):
    # TODO: Write Makefile to output directory, which, when executed, will perform species separation
    with fw.writing_to_file(fw.MakefileWriter, options[OUTPUT_DIR], "Makefile") as writer:
        _write_variable_definitions(logger, writer, options)
        _write_target_variable_definitions(logger, writer, options)
        _write_phony_targets(logger, writer, options)
        _write_all_target(logger, writer, options)
        _write_filtered_reads_target(logger, writer, options)
        _write_sorted_reads_target(logger, writer, options)
        _write_mapped_reads_target(logger, writer, options)
        _write_masked_reads_target(logger, writer, options)
        _write_collate_raw_reads_target(logger, writer, options)
        _write_mask_star_index_targets(logger, writer, options)
        _write_main_star_index_targets(logger, writer, options)
        _write_clean_target(logger, writer, options)


def _run_species_separation(logger, options):
    # TODO: Execute Makefile with nohup
    pass


def _separate_species(logger, options):
    os.mkdir(options[OUTPUT_DIR])

    _write_makefile(logger, options)

    if options[RUN_SEPARATION]:
        _run_species_separation(logger, options)


def separate_species(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="species_separator v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    _separate_species(logger, options)

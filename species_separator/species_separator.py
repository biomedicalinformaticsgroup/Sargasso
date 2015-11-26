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

e.g. bin/species_separator --reads-base-dir=/srv/data/ghardingham/neuron_astrocyte_activity/rnaseq --s1-index=/srv/data/genome/mouse/ensembl-80/STAR_Indices/primary_assembly --s2-index=/srv/data/genome/rat/ensembl-80/STAR_Indices/top_level -t 4 mouse rat samples.tsv my_results
"""

import docopt
import os
import os.path
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

SPECIES_NAME = "species-name"
GTF_FILE = "gtf-file"
GENOME_FASTA = "genome-fasta"
STAR_INDEX = "star-index"


# TODO: deal with single-end reads
class SampleInfo(object):
    """
    Encapsulates sample names and their accompanying reads files.
    """
    def __init__(self, base_reads_dir):
        """
        Create object.

        base_reads_dir: Common base path to all raw reads files.
        """
        self.base_reads_dir = base_reads_dir
        self.left_reads = {}
        self.right_reads = {}

    def get_sample_names(self):
        """
        Return a list of sample names.
        """
        return self.left_reads.keys()

    def get_left_reads(self, sample_name):
        """
        Return list of first in pair reads files for sample.
        """
        return self.left_reads[sample_name]

    def get_right_reads(self, sample_name):
        """
        Return list of second in pair reads files for sample.
        """
        return self.right_reads[sample_name]

    def add_sample_data(self, sample_data):
        """
        Add information for a single sample.

        sample_data: list of strings - sample name, comma-separated list of
        first in pair reads files, comma-separated list of second in pair reads
        files.
        """
        sample_name = sample_data[0]
        self.left_reads[sample_name] = sample_data[1].split(',')
        self.right_reads[sample_name] = sample_data[2].split(',')

    def validate(self):
        """
        Validate all raw reads files exist.
        """
        for reads_set in [self.left_reads, self.right_reads]:
            for reads_file_list in reads_set.values():
                for reads_file in reads_file_list:
                    self.validate_read_file(reads_file)

    def validate_read_file(self, reads_file):
        """
        Validate a particular raw reads file exists.

        reads_file: Path to a particular raw reads file (possibly truncated,
        and needing to be joined with 'base_reads_dir').
        """
        if self.base_reads_dir:
            reads_file = os.path.join(self.base_reads_dir, reads_file)
        opt.validate_file_option(
            reads_file, "Could not open reads file")


def _get_species_options(options, name_option, gtf_file_option,
                         genome_fasta_option, star_index_option):
    """
    Return a dictionary containing command-line options for a species.

    options: dictionary containing all command-line options.
    name_option: name of option specifying species name.
    gtf_file_option: name of option specifying GTF annotation file for the
    species.
    genome_fasta_option: name of option specifying directory containing genome
    FASTA files for the species.
    star_index_option: name of option specifying STAR index directory for the
    species.
    """
    species_options = {}
    species_options[SPECIES_NAME] = options[name_option]
    species_options[GTF_FILE] = options[gtf_file_option]
    species_options[GENOME_FASTA] = options[genome_fasta_option]
    species_options[STAR_INDEX] = options[star_index_option]
    return species_options


def _validate_species_options(species, species_options):
    """
    Validate command-line options for a species are correctly specified.

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


def _read_sample_info(options):
    """
    Return an object encapsulating samples and their accompanying read files.

    options: dictionary of command-line options
    """
    sample_info = SampleInfo(options[READS_BASE_DIR])

    for line in open(options[SAMPLES_FILE], 'r'):
        sample_data = line.split()

        if len(sample_data) != 3:
            line = line.rstrip('\n')
            if len(line) > 80:
                line = line[0:77] + "..."
            raise schema.SchemaError(
                None, "Sample file line should contain sample name and " +
                "lists of first and second pairs of read files, separated " +
                "by whitespace: \n{info}".format(info=line))

        sample_info.add_sample_data(sample_data)

    return sample_info


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
            options, SPECIES_ONE, SPECIES_ONE_GTF,
            SPECIES_ONE_GENOME_FASTA, SPECIES_ONE_INDEX)
        species_two_options = _get_species_options(
            options, SPECIES_TWO, SPECIES_TWO_GTF,
            SPECIES_TWO_GENOME_FASTA, SPECIES_TWO_INDEX)

        _validate_species_options("one", species_one_options)
        _validate_species_options("two", species_two_options)

        sample_info = _read_sample_info(options)
        sample_info.validate()

        return sample_info
    except schema.SchemaError as exc:
        # TODO: format exit message for 80 columns
        exit("Exiting: " + exc.code)


def _write_species_variable_definitions(
        logger, writer, species, species_options):
    """
    Write variable definitions for a particular species.

    logger: logging object
    writer: Makefile writer object
    species: species number string
    species_options: dictionary of options specific to a particular species.
    """
    writer.set_variable(
        "SPECIES_{species}".format(species=species),
        species_options[SPECIES_NAME])

    if species_options[GTF_FILE] is None:
        writer.set_variable(
            "SPECIES_{species}_GENOME_DIR".format(species=species),
            species_options[STAR_INDEX])
    else:
        writer.set_variable(
            "SPECIES_{species}_GTF".format(species=species),
            species_options[GTF_FILE])
        writer.set_variable(
            "SPECIES_{species}_GENOME_FASTA".format(species=species),
            species_options[GENOME_FASTA])

    writer.add_blank_line()


def _write_variable_definitions(logger, writer, options, sample_info):
    """
    Write variable definitions to Makefile.

    logger: logging object
    writer: Makefile writer object
    options: dictionary of command-line options
    sample_info: object encapsulating samples and their accompanying read files
    """
    writer.set_variable("NUM_THREADS", options[NUM_THREADS])
    writer.add_blank_line()

    sample_names = sample_info.get_sample_names()
    writer.set_variable("SAMPLES", " ".join(sample_names))
    writer.set_variable(
        "RAW_READS_DIRECTORY",
        options[READS_BASE_DIR] if options[READS_BASE_DIR] else "/")
    writer.set_variable(
        "RAW_READ_FILES_1",
        " ".join([",".join(sample_info.get_left_reads(name))
                  for name in sample_names]))
    writer.set_variable(
        "RAW_READ_FILES_2",
        " ".join([",".join(sample_info.get_right_reads(name))
                  for name in sample_names]))
    writer.add_blank_line()

    species_one_options = _get_species_options(
        options, SPECIES_ONE, SPECIES_ONE_GTF,
        SPECIES_ONE_GENOME_FASTA, SPECIES_ONE_INDEX)
    species_two_options = _get_species_options(
        options, SPECIES_TWO, SPECIES_TWO_GTF,
        SPECIES_TWO_GENOME_FASTA, SPECIES_TWO_INDEX)

    _write_species_variable_definitions(
        logger, writer, "ONE", species_one_options)
    _write_species_variable_definitions(
        logger, writer, "TWO", species_two_options)


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


def _write_makefile(logger, options, sample_info):
    # TODO: Write Makefile to output directory, which, when executed, will perform species separation
    with fw.writing_to_file(fw.MakefileWriter, options[OUTPUT_DIR], "Makefile") as writer:
        _write_variable_definitions(logger, writer, options, sample_info)
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


def _separate_species(logger, options, sample_info):
    os.mkdir(options[OUTPUT_DIR])

    _write_makefile(logger, options, sample_info)

    if options[RUN_SEPARATION]:
        _run_species_separation(logger, options)


def separate_species(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="species_separator v" + __version__)

    # Validate command-line options
    sample_info = _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    _separate_species(logger, options, sample_info)

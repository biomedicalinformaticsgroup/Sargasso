#!/usr/bin/env python

"""Usage:
    species_separator [--log-level=<log-level>] [--reads-base-dir=<reads-base-dir>] [--num-threads=<num-threads>] [--s1-gtf=<species-one-gtf-file>] [--s2-gtf=<species-two-gtf-file>] [--s1-genome-fasta=<species-one-genome-fasta>] [--s2-genome-fasta=<species-two-genome-fasta>] [--s1-index=<species-one-star-index>] [--s2-index=<species-two-star-index>] [--mismatch-threshold=<mismatch-threshold>] --minmatch-threshold=<minmatch-threshold> [--multimap-threshold=<multimap-threshold>] [--run-separation] <species-one> <species-two> <samples-file> <output-dir>

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
--mismatch-threshold=<mismatch-threshold>       Maximum number of mismatches allowed during filtering [default: 0].
--minmatch-threshold=<minmatch-threshold>       Minimum number of read bases that must be perfectly matched
--multimap-threshold=<multimap-threshold>       Maximum number of multiple mappings allowed during filtering [default: 1].
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
from . import process
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
MISMATCH_THRESHOLD = "--mismatch-threshold"
MINMATCH_THRESHOLD = "--minmatch-threshold"
MULTIMAP_THRESHOLD = "--multimap-threshold"
RUN_SEPARATION = "--run-separation"

SPECIES_NAME = "species-name"
GTF_FILE = "gtf-file"
GENOME_FASTA = "genome-fasta"
STAR_INDEX = "star-index"

ALL_TARGET = "all"
CLEAN_TARGET = "clean"
STAR_INDICES_TARGET = "STAR_INDICES"
COLLATE_RAW_READS_TARGET = "COLLATE_RAW_READS"
MAPPED_READS_TARGET = "MAPPED_READS"
SORTED_READS_TARGET = "SORTED_READS"
FILTERED_READS_TARGET = "FILTERED_READS"

NUM_THREADS_VARIABLE = "NUM_THREADS"
SPECIES_ONE_VARIABLE = "SPECIES_ONE"
SPECIES_TWO_VARIABLE = "SPECIES_TWO"
SPECIES_GTF_FILE_VARIABLE = "{species}_GTF"
SPECIES_GENOME_FASTA_VARIABLE = "{species}_GENOME_FASTA"
SPECIES_STAR_INDEX_VARIABLE = "{species}_GENOME_DIR"
SAMPLES_VARIABLE = "SAMPLES"
RAW_READS_DIRECTORY_VARIABLE = "RAW_READS_DIRECTORY"
RAW_READS_LEFT_VARIABLE = "RAW_READS_FILES_1"
RAW_READS_RIGHT_VARIABLE = "RAW_READS_FILES_2"

SINGLE_END_READS_TYPE = "single"
PAIRED_END_READS_TYPE = "paired"


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
        self.paired_end = None

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

    def paired_end_reads(self):
        """
        Returns True iff sample read files are paired-end data.
        """
        return self.paired_end

    def add_sample_data(self, sample_data):
        """
        Add information for a single sample.

        sample_data: list of strings - sample name, comma-separated list of
        read files (or comma-separated list of first in pair reads files,
        comma-separated list of second in pair reads files).
        """
        sample_name = sample_data[0]
        self.left_reads[sample_name] = sample_data[1].split(',')

        if len(sample_data) == 3:
            self.right_reads[sample_name] = sample_data[2].split(',')
            self.paired_end = True
        else:
            self.paired_end = False

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

        if len(sample_data) < 2 or len(sample_data) > 3:
            line = line.rstrip('\n')
            if len(line) > 80:
                line = line[0:77] + "..."
            raise schema.SchemaError(
                None, "Sample file line should contain sample name and " +
                "lists of read files (or first and second pairs of read " +
                "files), separated by whitespace: \n{info}".format(info=line))

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
        options[MISMATCH_THRESHOLD] = opt.validate_int_option(
            options[MISMATCH_THRESHOLD],
            "Maximum number of mismatches must be a non-negative integer.",
            0, True)
        options[MINMATCH_THRESHOLD] = opt.validate_int_option(
            options[MINMATCH_THRESHOLD],
            "Minimum number of perfect matches must be a positive integer.",
            1)
        options[MULTIMAP_THRESHOLD] = opt.validate_int_option(
            options[MULTIMAP_THRESHOLD],
            "Maximum number of multiple mappings must be a positive integer.",
            1)

        species_one_options = _get_species_options(
            options, SPECIES_ONE, SPECIES_ONE_GTF,
            SPECIES_ONE_GENOME_FASTA, SPECIES_ONE_INDEX)
        species_two_options = _get_species_options(
            options, SPECIES_TWO, SPECIES_TWO_GTF,
            SPECIES_TWO_GENOME_FASTA, SPECIES_TWO_INDEX)

        _validate_species_options("one", species_one_options)
        _validate_species_options("two", species_two_options)

        sample_info = _read_sample_info(options)
        # TODO: validate that all samples consistently have either single- or
        # paired-end reads
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
        "{species}".format(species=species),
        species_options[SPECIES_NAME])

    if species_options[GTF_FILE] is None:
        writer.set_variable(
            "{species}_GENOME_DIR".format(species=species),
            species_options[STAR_INDEX])
    else:
        writer.set_variable(
            "{species}_GTF".format(species=species),
            species_options[GTF_FILE])
        writer.set_variable(
            "{species}_GENOME_FASTA".format(species=species),
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
    writer.set_variable(NUM_THREADS_VARIABLE, options[NUM_THREADS])
    writer.add_blank_line()

    sample_names = sample_info.get_sample_names()
    writer.set_variable(SAMPLES_VARIABLE, " ".join(sample_names))
    writer.set_variable(
        RAW_READS_DIRECTORY_VARIABLE,
        options[READS_BASE_DIR] if options[READS_BASE_DIR] else "/")
    writer.set_variable(
        RAW_READS_LEFT_VARIABLE,
        " ".join([",".join(sample_info.get_left_reads(name))
                  for name in sample_names]))

    if sample_info.paired_end_reads():
        writer.set_variable(
            RAW_READS_RIGHT_VARIABLE,
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
        logger, writer, SPECIES_ONE_VARIABLE, species_one_options)
    _write_species_variable_definitions(
        logger, writer, SPECIES_TWO_VARIABLE, species_two_options)


def _write_target_variable_definitions(logger, writer):
    """
    Write target directory variable definitions to Makefile.

    logger: logging object
    writer: Makefile writer object
    """
    writer.set_variable(STAR_INDICES_TARGET, "star_indices")
    writer.set_variable(COLLATE_RAW_READS_TARGET, "raw_reads")
    writer.set_variable(MAPPED_READS_TARGET, "mapped_reads")
    writer.set_variable(SORTED_READS_TARGET, "sorted_reads")
    writer.set_variable(FILTERED_READS_TARGET, "filtered_reads")
    writer.add_blank_line()


def _write_phony_targets(logger, writer):
    """
    Write phony target definitions to Makefile.

    logger: logging object
    writer: Makefile writer object
    """
    with writer.target_definition(
            ".PHONY", [ALL_TARGET, CLEAN_TARGET],
            raw_target=True, raw_dependencies=True):
        pass


def _write_all_target(logger, writer):
    """
    Write main target definition to Makefile.

    logger: logging object
    writer: Makefile writer object
    """
    with writer.target_definition(
            ALL_TARGET, [FILTERED_READS_TARGET], raw_target=True):
        pass


def _write_filtered_reads_target(logger, writer, options):
    """
    Write target to separate reads by species to Makefile.

    logger: logging object
    writer: Makefile writer object
    """
    with writer.target_definition(
            FILTERED_READS_TARGET, [SORTED_READS_TARGET]):
        writer.add_comment(
            "For each sample, take the reads mapping to each genome and " +
            "filter them to their correct species of origin")
        writer.make_target_directory(FILTERED_READS_TARGET)
        writer.add_command(
            "filter_reads",
            [writer.variable_val(SPECIES_ONE_VARIABLE),
             writer.variable_val(SPECIES_TWO_VARIABLE),
             "\"{var}\"".format(var=writer.variable_val(SAMPLES_VARIABLE)),
             writer.variable_val(SORTED_READS_TARGET),
             writer.variable_val(FILTERED_READS_TARGET),
             writer.variable_val(NUM_THREADS_VARIABLE),
             options[MISMATCH_THRESHOLD],
             options[MINMATCH_THRESHOLD],
             options[MULTIMAP_THRESHOLD]])


def _write_sorted_reads_target(logger, writer):
    """
    Write target to sort reads by name to Makefile.

    logger: logging object
    writer: Makefile writer object
    """
    with writer.target_definition(
            SORTED_READS_TARGET, [MAPPED_READS_TARGET]):
        writer.add_comment(
            "For each sample, sort the mapped reads into read name order")
        writer.make_target_directory(SORTED_READS_TARGET)
        writer.add_command(
            "sort_reads",
            ["\"{s1} {s2}\"".format(
                s1=writer.variable_val(SPECIES_ONE_VARIABLE),
                s2=writer.variable_val(SPECIES_TWO_VARIABLE)),
             "\"{var}\"".format(
                 var=writer.variable_val(SAMPLES_VARIABLE)),
             writer.variable_val(NUM_THREADS_VARIABLE),
             writer.variable_val(MAPPED_READS_TARGET),
             writer.variable_val(SORTED_READS_TARGET)])


def _write_mapped_reads_target(logger, writer, sample_info):
    """
    Write target to map reads to each species to Makefile.

    logger: logging object
    writer: Makefile writer object
    """
    with writer.target_definition(
        MAPPED_READS_TARGET,
        ["{index}/{s1}".format(index=writer.variable_val(STAR_INDICES_TARGET),
                               s1=writer.variable_val(SPECIES_ONE_VARIABLE)),
         "{index}/{s2}".format(index=writer.variable_val(STAR_INDICES_TARGET),
                               s2=writer.variable_val(SPECIES_TWO_VARIABLE)),
         writer.variable_val(COLLATE_RAW_READS_TARGET)],
            raw_dependencies=True):
        writer.add_comment(
            "Map reads for each sample to each species' genome")
        writer.make_target_directory(MAPPED_READS_TARGET)

        map_reads_params = [
            "\"{s1} {s2}\"".format(
                s1=writer.variable_val(SPECIES_ONE_VARIABLE),
                s2=writer.variable_val(SPECIES_TWO_VARIABLE)),
            "\"{var}\"".format(var=writer.variable_val(SAMPLES_VARIABLE)),
            writer.variable_val(STAR_INDICES_TARGET),
            writer.variable_val(NUM_THREADS_VARIABLE),
            writer.variable_val(COLLATE_RAW_READS_TARGET),
            writer.variable_val(MAPPED_READS_TARGET)]

        map_reads_params.append(
            PAIRED_END_READS_TYPE if sample_info.paired_end_reads()
            else SINGLE_END_READS_TYPE)

        writer.add_command("map_reads", map_reads_params)


#def _write_masked_reads_target(logger, writer, options):
    # TODO: Write "masked_reads" target
    #pass


def _write_collate_raw_reads_target(logger, writer, sample_info):
    """
    Write target to collect raw reads files to Makefile.

    logger: logging object
    writer: Makefile writer object
    sample_info: object encapsulating samples and their accompanying read files
    """
    with writer.target_definition(COLLATE_RAW_READS_TARGET, []):
        writer.add_comment(
            "Create a directory with sub-directories for each sample, " +
            "each of which contains links to the input raw reads files " +
            "for that sample")
        writer.make_target_directory(COLLATE_RAW_READS_TARGET)

        collate_raw_reads_params = [
            "\"{var}\"".format(var=writer.variable_val(SAMPLES_VARIABLE)),
            writer.variable_val(RAW_READS_DIRECTORY_VARIABLE),
            writer.variable_val(COLLATE_RAW_READS_TARGET),
        ]

        if sample_info.paired_end_reads():
            collate_raw_reads_params += [
                PAIRED_END_READS_TYPE,
                "\"{var}\"".format(
                    var=writer.variable_val(RAW_READS_LEFT_VARIABLE)),
                "\"{var}\"".format(
                    var=writer.variable_val(RAW_READS_RIGHT_VARIABLE))
            ]
        else:
            collate_raw_reads_params += [
                SINGLE_END_READS_TYPE,
                "\"{var}\"".format(
                    var=writer.variable_val(RAW_READS_LEFT_VARIABLE)),
                "\"\""
            ]

        writer.add_command("collate_raw_reads", collate_raw_reads_params)


#def _write_mask_star_index_targets(logger, writer, options):
    # TODO: Write targets for STAR indices for mask sequences
    #pass


def _write_species_main_star_index_target(
        logger, writer, species_var, species_options):
    """
    Write target to create or link to STAR index for a species to Makefile.

    logger: logging object
    writer: Makefile writer object
    species_var: Makefile variable for species
    species_options: dictionary of command-line options for the particular
    species
    """
    target = "{index}/{spec}".format(
        index=writer.variable_val(STAR_INDICES_TARGET),
        spec=writer.variable_val(species_var))

    with writer.target_definition(target, [], raw_target=True):
        writer.make_target_directory(STAR_INDICES_TARGET)

        if species_options[GTF_FILE] is None:
            writer.add_command(
                "ln",
                ["-s",
                 writer.variable_val(
                     SPECIES_STAR_INDEX_VARIABLE.format(species=species_var)),
                 target])
        else:
            writer.add_command(
                "build_star_index",
                [writer.variable_val(
                    SPECIES_GENOME_FASTA_VARIABLE.format(species=species_var)),
                 writer.variable_val(
                     SPECIES_GTF_FILE_VARIABLE.format(species=species_var)),
                 writer.variable_val(NUM_THREADS_VARIABLE),
                 target])


def _write_main_star_index_targets(logger, writer, options):
    """
    Write targets to create or link to STAR indices to Makefile.

    logger: logging object
    writer: Makefile writer object
    options: dictionary of command-line options
    """
    species_one_options = _get_species_options(
        options, SPECIES_ONE, SPECIES_ONE_GTF,
        SPECIES_ONE_GENOME_FASTA, SPECIES_ONE_INDEX)
    species_two_options = _get_species_options(
        options, SPECIES_TWO, SPECIES_TWO_GTF,
        SPECIES_TWO_GENOME_FASTA, SPECIES_TWO_INDEX)

    _write_species_main_star_index_target(
        logger, writer, SPECIES_ONE_VARIABLE, species_one_options)
    _write_species_main_star_index_target(
        logger, writer, SPECIES_TWO_VARIABLE, species_two_options)


def _write_clean_target(logger, writer):
    """
    Write target to clean results directory to Makefile.

    logger: logging object
    writer: Makefile writer object
    """
    with writer.target_definition(CLEAN_TARGET, [], raw_target=True):
        writer.remove_target_directory(STAR_INDICES_TARGET)
        writer.remove_target_directory(COLLATE_RAW_READS_TARGET)
        writer.remove_target_directory(MAPPED_READS_TARGET)
        writer.remove_target_directory(SORTED_READS_TARGET)
        writer.remove_target_directory(FILTERED_READS_TARGET)


def _write_makefile(logger, options, sample_info):
    """
    Write Makefile to results directory to perform species separation.

    logger: logging object
    options: dictionary of command-line options
    sample_info: object encapsulating samples and their accompanying read files
    """
    with fw.writing_to_file(
            fw.MakefileWriter, options[OUTPUT_DIR], "Makefile") as writer:
        _write_variable_definitions(logger, writer, options, sample_info)
        _write_target_variable_definitions(logger, writer)
        _write_phony_targets(logger, writer)
        _write_all_target(logger, writer)
        _write_filtered_reads_target(logger, writer, options)
        _write_sorted_reads_target(logger, writer)
        _write_mapped_reads_target(logger, writer, sample_info)
        #_write_masked_reads_target(logger, writer)
        _write_collate_raw_reads_target(logger, writer, sample_info)
        #_write_mask_star_index_targets(logger, writer, options)
        _write_main_star_index_targets(logger, writer, options)
        _write_clean_target(logger, writer)


def _run_species_separation(logger, options):
    """
    Executes the written Makefile with nohup.

    logger: logging object
    options: dictionary of command-line options
    """
    process.run_in_directory(options[OUTPUT_DIR], "make")


def separate_species(args):
    """
    Write a Makefile to perform species separation and optionally execute it.

    args: list of command-line arguments
    """
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="species_separator v" + __version__)

    # Validate command-line options
    sample_info = _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Create output directory
    os.mkdir(options[OUTPUT_DIR])

    # Write Makefile to output directory
    _write_makefile(logger, options, sample_info)

    # If specified, execute the Makefile with nohup
    if options[RUN_SEPARATION]:
        _run_species_separation(logger, options)

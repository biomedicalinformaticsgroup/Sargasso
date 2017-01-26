#!/usr/bin/env python

"""Usage:
    species_separator
        [--log-level=<log-level>]
        [--reads-base-dir=<reads-base-dir>] [--num-threads=<num-threads>]
        [--mismatch-threshold=<mismatch-threshold>]
        [--minmatch-threshold=<minmatch-threshold>]
        [--multimap-threshold=<multimap-threshold>]
        [--overhang-threshold=<overhang-threshold>]
        [--reject-multimaps] [--best] [--conservative] [--recall]
        [--run-separation]
        <samples-file> <output-dir>
        (<species> <species-star-info>)
        (<species> <species-star-info>)
        ...

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<samples-file>
    TSV file giving paths (relative to <reads-base-dir>) of raw RNA-seq read
    data files for each sample.
<output-dir>
    Output directory into which Makefile will be written, and in which species
    separation will be performed.
<species>
    Name of species.
<species-star-info>
    Either a STAR index directory for the species, or a comma-separated list of
    (i) a GTF annotation file and (ii) a directory containing genome FASTA
    files for the species.
<species-genome-fasta>
    Directory containing genome FASTA files for species.
<species-index>
    STAR index directory for species.
--reads-base-dir=<reads-base-dir>
    Base directory for raw RNA-seq read data files.
-t <num-threads> --num-threads=<num-threads>
    Number of threads to use for parallel processing. [default: 1]
--mismatch-threshold=<mismatch-threshold>
    Maximum percentage of bases allowed to be mismatches against the genome
    during filtering. For single-end reads, the total number of bases is the
    read length; for paired-end reads it is twice the read length [default: 0].
--minmatch-threshold=<minmatch-threshold>
    Maximum percentage of read length allowed to not be mapped during
    filtering. Read length refers to the length of a single read in both
    single- and paired-end cases [default: 0].
--multimap-threshold=<multimap-threshold>
    Maximum number of multiple mappings allowed during filtering [default: 1].
--overhang-threshold=<overhang-threshold>
    If set, allows specification of the minimum number of bases that are
    allowed on either side of an exon boundary for a read mapping to be
    accepted [default: 5].
--reject-multimaps
    If set, any read which multimaps to either species' genome will be rejected
    and not be assigned to either species.
--best
    Adopt a filtering strategy that provides an excellent balance between
    sensitivity and specificity. Note that specifying this option overrides the
    values of the mismatch-threshold, minmatch-threshold and
    multimap-threshold options. In addition, reject-multimaps is turned off.
--conservative
    Adopt a filtering strategy where minimising the number of reads
    mis-assigned to the wrong species takes foremost priority. Note that
    specifying this option overrides the values of the mismatch-threshold,
    minmatch-threshold and multimap-threshold options. In addition,
    reject-multimaps is turned on.
--recall
    Adopt a filtering strategy where sensitivity is prioritised over
    specificity.  Note that specifying this option overrides the values of the
    mismatch-threshold, minmatch-threshold and multimap-threshold options. In
    addition, reject-multimaps is turned off.
--run-separation
    If specified, species separation will be run; otherwise scripts to perform
    separation will be created but not run.

Given a set of RNA-seq samples containing mixed-species read data, determine,
where possible, from which of the species each read originated. Mapped
reads are written to per-sample and -species specific output BAM files.

Species separation of mixed-species RNA-seq data is performed in a number of
stages, of which the most important steps are:

1) Mapping of raw RNA-seq data to each species' genome using the STAR read
aligner.
2) Sorting of mapped RNA-seq data.
3) Assignment of mapped, sorted RNA-seq reads to their correct species of
origin.

If the option "--run-separation" is not specified, a Makefile is written to the
given output directory, via which all stages of species separation can be
run under the user's control. If "--run-separation" is specified however, the
Makefile is both written and executed, and all stages of species separation are
performed automatically.

n.b. Many stages of species separation can be executed across multiple threads
by specifying the "--num-threads" option.

e.g.:

species_separator --reads-base-dir=/srv/data/rnaseq --num-threads 4 --run-separation samples.tsv my_results mouse /srv/data/genome/mouse/STAR_Index rat /srv/data/genome/rat/STAR_Index
"""

import docopt
import os
import os.path
import schema

from . import file_writer as fw
from . import filter_sample_reads
from . import options as opt
from . import process
from .__init__ import __version__
from datetime import datetime

SAMPLES_FILE = "<samples-file>"
OUTPUT_DIR = "<output-dir>"
SPECIES = "<species>"
SPECIES_STAR_INFO = "<species-star-info>"
READS_BASE_DIR = "--reads-base-dir"
NUM_THREADS = "--num-threads"
MISMATCH_THRESHOLD = "--mismatch-threshold"
MINMATCH_THRESHOLD = "--minmatch-threshold"
MULTIMAP_THRESHOLD = "--multimap-threshold"
OVERHANG_THRESHOLD = "--overhang-threshold"
REJECT_MULTIMAPS = "--reject-multimaps"
OPTIMAL_STRATEGY = "--best"
CONSERVATIVE_STRATEGY = "--conservative"
RECALL_STRATEGY = "--recall"
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
SAMPLES_VARIABLE = "SAMPLES"
RAW_READS_DIRECTORY_VARIABLE = "RAW_READS_DIRECTORY"
RAW_READS_LEFT_VARIABLE = "RAW_READS_FILES_1"
RAW_READS_RIGHT_VARIABLE = "RAW_READS_FILES_2"

SINGLE_END_READS_TYPE = "single"
PAIRED_END_READS_TYPE = "paired"

EXECUTION_RECORD_ENTRIES = [
    ["Samples File", SAMPLES_FILE],
    ["Output Dir", OUTPUT_DIR],
    ["Species", SPECIES],
    ["Species STAR info", SPECIES_STAR_INFO],
    ["Reads Base Dir", READS_BASE_DIR],
    ["Number of Threads", NUM_THREADS],
    ["Mismatch Threshold", MISMATCH_THRESHOLD],
    ["Minmatch Threshold", MINMATCH_THRESHOLD],
    ["Multimap Threshold", MULTIMAP_THRESHOLD],
    ["Overhang Threshold", OVERHANG_THRESHOLD],
    ["Reject Multimaps", REJECT_MULTIMAPS],
    ["Optimal Strategy", OPTIMAL_STRATEGY],
    ["Conservative Strategy", CONSERVATIVE_STRATEGY],
    ["Recall Strategy", RECALL_STRATEGY],
    ["Run Separation", RUN_SEPARATION],
]


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


def _get_species_options(options, species_index):
    """
    Return a dictionary containing command-line options for a species.

    options: dictionary containing all command-line options.
    species_index: which species to return options for
    """

    species_options = { SPECIES_NAME: options[SPECIES][species_index] }

    star_infos = options[SPECIES_STAR_INFO][species_index].split(",")

    if len(star_infos) == 1:
        species_options[STAR_INDEX] = star_infos[0]
        species_options[GTF_FILE] = None
        species_options[GENOME_FASTA] = None
    elif len(star_infos) == 2:
        species_options[STAR_INDEX] = None
        species_options[GTF_FILE] = star_infos[0]
        species_options[GENOME_FASTA] = star_infos[1]
    else:
        raise schema.SchemaError(
            None, "Should specify either a STAR index or both GTF file " +
                  "and genome FASTA directory for species {species}.".
                    format(species=species_options[SPECIES_NAME]))

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

        if options[OPTIMAL_STRATEGY]:
            options[MISMATCH_THRESHOLD] = 1
            options[MINMATCH_THRESHOLD] = 2
            options[MULTIMAP_THRESHOLD] = 999999
            options[REJECT_MULTIMAPS] = False
        elif options[CONSERVATIVE_STRATEGY]:
            options[MISMATCH_THRESHOLD] = 0
            options[MINMATCH_THRESHOLD] = 0
            options[MULTIMAP_THRESHOLD] = 1
            options[REJECT_MULTIMAPS] = True
        elif options[RECALL_STRATEGY]:
            options[MISMATCH_THRESHOLD] = 2
            options[MINMATCH_THRESHOLD] = 10
            options[MULTIMAP_THRESHOLD] = 999999
            options[REJECT_MULTIMAPS] = False

        filter_sample_reads.validate_threshold_options(
            options, MISMATCH_THRESHOLD, MINMATCH_THRESHOLD,
            MULTIMAP_THRESHOLD, OVERHANG_THRESHOLD)

        for i, species in enumerate(options[SPECIES]):
            species_options = _get_species_options(options, i)
            _validate_species_options(species, species_options)

        sample_info = _read_sample_info(options)
        # TODO: validate that all samples consistently have either single- or
        # paired-end reads
        sample_info.validate()

        return sample_info
    except schema.SchemaError as exc:
        exit("Exiting: " + exc.code)


def _get_gtf_file_variable(species):
    """
    Return GTF file variable name for a species.

    species: species identifier
    """
    return "{species}_GTF".format(species=species.upper())


def _get_genome_fasta_variable(species):
    """
    Return genome FASTA variable name for a species.

    species: species identifier
    """
    return "{species}_GENOME_FASTA_DIR".format(species=species.upper())


def _get_star_index_variable(species):
    """
    Return STAR index variable name for a species.

    species: species identifier
    """
    return "{species}_STAR_DIR".format(species=species.upper())


def _write_species_variable_definitions(
        logger, writer, species, species_options):
    """
    Write variable definitions for a particular species.

    logger: logging object
    writer: Makefile writer object
    species: species number string
    species_options: dictionary of options specific to a particular species.
    """
    if species_options[GTF_FILE] is None:
        writer.set_variable(
            _get_star_index_variable(species), species_options[STAR_INDEX])
    else:
        writer.set_variable(
            _get_gtf_file_variable(species), species_options[GTF_FILE])
        writer.set_variable(
            _get_genome_fasta_variable(species), species_options[GENOME_FASTA])


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

    for i, species in enumerate(options[SPECIES]):
        species_options = _get_species_options(options, i)
        _write_species_variable_definitions(
            logger, writer, species, species_options)

    writer.add_blank_line()


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

        writer.add_command("filter_reads", [
             "\"{var}\"".format(var=writer.variable_val(SAMPLES_VARIABLE)),
             writer.variable_val(SORTED_READS_TARGET),
             writer.variable_val(FILTERED_READS_TARGET),
             writer.variable_val(NUM_THREADS_VARIABLE),
             options[MISMATCH_THRESHOLD],
             options[MINMATCH_THRESHOLD],
             options[MULTIMAP_THRESHOLD],
             options[OVERHANG_THRESHOLD],
             "--reject-multimaps" if options[REJECT_MULTIMAPS] else "\"\"",
             "{sl}".format(sl=" ".join(options[SPECIES]))])


def _write_sorted_reads_target(logger, writer, options):
    """
    Write target to sort reads by name to Makefile.

    logger: logging object
    writer: Makefile writer object
    options: dictionary of command-line options
    """
    with writer.target_definition(
            SORTED_READS_TARGET, [MAPPED_READS_TARGET]):
        writer.add_comment(
            "For each sample, sort the mapped reads into read name order")
        writer.make_target_directory(SORTED_READS_TARGET)

        writer.add_command(
            "sort_reads",
            ["\"{sl}\"".format(sl=" ".join(options[SPECIES])),
             "\"{var}\"".format(
                 var=writer.variable_val(SAMPLES_VARIABLE)),
             writer.variable_val(NUM_THREADS_VARIABLE),
             writer.variable_val(MAPPED_READS_TARGET),
             writer.variable_val(SORTED_READS_TARGET)])


def _write_mapped_reads_target(logger, writer, sample_info, options):
    """
    Write target to map reads to each species to Makefile.

    logger: logging object
    writer: Makefile writer object
    """

    index_targets = ["{index}/{species}".format(
                        index=writer.variable_val(STAR_INDICES_TARGET),
                        species=species)
                     for species in options[SPECIES]]

    index_targets.append(writer.variable_val(COLLATE_RAW_READS_TARGET))

    with writer.target_definition(
            MAPPED_READS_TARGET, index_targets, raw_dependencies=True):
        writer.add_comment(
            "Map reads for each sample to each species' genome")
        writer.make_target_directory(MAPPED_READS_TARGET)

        map_reads_params = [
            "\"{sl}\"".format(sl=" ".join(options[SPECIES])),
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
        logger, writer, species, species_options):
    """
    Write target to create or link to STAR index for a species to Makefile.

    logger: logging object
    writer: Makefile writer object
    species_var: Makefile variable for species
    species_options: dictionary of command-line options for the particular
    species
    """
    target = "{index}/{spec}".format(
        index=writer.variable_val(STAR_INDICES_TARGET), spec=species)

    with writer.target_definition(target, [], raw_target=True):
        if species_options[GTF_FILE] is None:
            writer.make_target_directory(STAR_INDICES_TARGET)
            writer.add_command(
                "ln",
                ["-s",
                 writer.variable_val(_get_star_index_variable(species)),
                 target])
        else:
            writer.make_target_directory(target, raw_target=True)
            writer.add_command(
                "build_star_index",
                [writer.variable_val(_get_genome_fasta_variable(species)),
                 writer.variable_val(_get_gtf_file_variable(species)),
                 writer.variable_val(NUM_THREADS_VARIABLE),
                 target])


def _write_main_star_index_targets(logger, writer, options):
    """
    Write targets to create or link to STAR indices to Makefile.

    logger: logging object
    writer: Makefile writer object
    options: dictionary of command-line options
    """
    for i, species in enumerate(options[SPECIES]):
        species_options = _get_species_options(options, i)
        _write_species_main_star_index_target(
            logger, writer, species, species_options)


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
        _write_sorted_reads_target(logger, writer, options)
        _write_mapped_reads_target(logger, writer, sample_info, options)
        #_write_masked_reads_target(logger, writer)
        _write_collate_raw_reads_target(logger, writer, sample_info)
        #_write_mask_star_index_targets(logger, writer, options)
        _write_main_star_index_targets(logger, writer, options)
        _write_clean_target(logger, writer)


def _write_execution_record(options):
    """
    Write a log file containing all execution parameters in addition to the
    time and date of execution

    options: dictionary of command-line options
    """
    out_text = "Execution Record - {t}\n".format(
        t=str(datetime.now().isoformat()))
    out_text += "\n".join(["{desc}: {val}".format(
                           desc=it[0], val=str(options[it[1]]))
                           for it in EXECUTION_RECORD_ENTRIES])

    out_file = os.path.join(options[OUTPUT_DIR], "execution_record.txt")
    with open(out_file, 'w') as erf:
        erf.write(out_text)


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

    # Write Execution Record to output directory
    _write_execution_record(options)

    # If specified, execute the Makefile with nohup
    if options[RUN_SEPARATION]:
        _run_species_separation(logger, options)

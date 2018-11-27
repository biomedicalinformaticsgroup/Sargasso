import logging


class Options(object):
    DATA_TYPE = "<data-type>"
    SAMPLES_FILE = "<samples-file>"
    OUTPUT_DIR = "<output-dir>"
    SPECIES = "<species>"
    SPECIES_INFO = "<species-info>"
    READS_BASE_DIR = "--reads-base-dir"
    NUM_THREADS = "--num-threads"
    MISMATCH_THRESHOLD = "--mismatch-threshold"
    MINMATCH_THRESHOLD = "--minmatch-threshold"
    MULTIMAP_THRESHOLD = "--multimap-threshold"
    REJECT_MULTIMAPS = "--reject-multimaps"
    OPTIMAL_STRATEGY = "--best"
    CONSERVATIVE_STRATEGY = "--conservative"
    RECALL_STRATEGY = "--recall"
    PERMISSIVE_STRATEGY = "--permissive"
    RUN_SEPARATION = "--run-separation"
    DELETE_INTERMEDIATE = "--delete-intermediate"
    MAPPER_EXECUTABLE = "--mapper-executable"
    MAPPER_INDEX_EXECUTABLE = "--mapper-index-executable"
    SAMBAMBA_SORT_TMP_DIR = "--sambamba-sort-tmp-dir"

    SPECIES_NAME = "species-name"
    GTF_FILE = "gtf-file"
    GENOME_FASTA = "genome-fasta"
    MAPPER_INDEX = "mapper-index"

    ALL_TARGET = "all"
    CLEAN_TARGET = "clean"
    MAPPER_INDICES_TARGET = "MAPPER_INDICES"
    COLLATE_RAW_READS_TARGET = "COLLATE_RAW_READS"
    MAPPED_READS_TARGET = "MAPPED_READS"
    SORTED_READS_TARGET = "SORTED_READS"
    FILTERED_READS_TARGET = "FILTERED_READS"

    DATA_TYPE_VARIABLE = "DATA_TYPE"
    NUM_THREADS_VARIABLE = "NUM_THREADS"
    MAPPER_EXECUTABLE_VARIABLE = "MAPPER_EXECUTABLE"
    SAMBAMBA_SORT_TMP_DIR_VARIABLE = "SAMBAMBA_SORT_TMP_DIR"
    SAMPLES_VARIABLE = "SAMPLES"
    RAW_READS_DIRECTORY_VARIABLE = "RAW_READS_DIRECTORY"
    RAW_READS_LEFT_VARIABLE = "RAW_READS_FILES_1"
    RAW_READS_RIGHT_VARIABLE = "RAW_READS_FILES_2"

    SINGLE_END_READS_TYPE = "single"
    PAIRED_END_READS_TYPE = "paired"

    SUPPORTED_DATATYPE = {'rnaseq': True, 'chipseq': True}

    SAMPLE_INFO_INDEX = "sample_info"
    SPECIES_OPTIONS_INDEX = "species_options"

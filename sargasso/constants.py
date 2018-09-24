#!/usr/bin/env python

SAMPLES_FILE = "<samples-file>"
OUTPUT_DIR = "<output-dir>"
SPECIES = "<species>"
SPECIES_STAR_INFO = "<species-star-info>"
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
STAR_EXECUTABLE = "--star-executable"
SAMBAMBA_SORT_TMP_DIR = "--sambamba-sort-tmp-dir"

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
STAR_EXECUTABLE_VARIABLE = "STAR_EXECUTABLE"
SAMBAMBA_SORT_TMP_DIR_VARIABLE = "SAMBAMBA_SORT_TMP_DIR"
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
    ["Reject Multimaps", REJECT_MULTIMAPS],
    ["Optimal Strategy", OPTIMAL_STRATEGY],
    ["Conservative Strategy", CONSERVATIVE_STRATEGY],
    ["Recall Strategy", RECALL_STRATEGY],
    ["Run Separation", RUN_SEPARATION],
    ["Delete Intermediate" , DELETE_INTERMEDIATE],
    ["Star Executable Path", STAR_EXECUTABLE],
    ["Sambamba Sort Tmp Dir", SAMBAMBA_SORT_TMP_DIR],
]



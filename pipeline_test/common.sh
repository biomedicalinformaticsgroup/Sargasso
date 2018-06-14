function test_expected_results {
    SPECIES_ID=$1
    EXPECTED_RESULTS_STRING=$2
    MISMATCH_THRESHOLD=$3
    MINMATCH_THRESHOLD=$4
    MULTIMAP_THRESHOLD=$5

    if ! grep --quiet "${EXPECTED_RESULTS_STRING}" ${LOG_FILE}; then
      echo "For mismatch-threshold=${MISMATCH_THRESHOLD}, minmatch-threshold=${MINMATCH_THRESHOLD}, multimap-threshold=${MULTIMAP_THRESHOLD}"
      echo "Expected: ${EXPECTED_RESULTS_STRING}"

      grep_res=$(grep "Species ${SPECIES_ID}:" ${LOG_FILE})
      echo "Got: ${grep_res}"
      exit 1
    fi
}

function expected_results {
    SPECIES_ID=$1
    FILTERED_HITS_EXPECTED=$2
    FILTERED_READS_EXPECTED=$3
    REJECTED_HITS_EXPECTED=$4
    REJECTED_READS_EXPECTED=$5
    AMBIGUOUS_HITS_EXPECTED=$6
    AMBIGUOUS_READS_EXPECTED=$7

    echo "Species ${SPECIES_ID}: wrote ${FILTERED_HITS_EXPECTED} filtered hits for ${FILTERED_READS_EXPECTED} reads; ${REJECTED_HITS_EXPECTED} hits for ${REJECTED_READS_EXPECTED} reads were rejected outright, and ${AMBIGUOUS_HITS_EXPECTED} hits for ${AMBIGUOUS_READS_EXPECTED} reads were rejected as ambiguous."
}

MAIN_DIR=$(pwd)

DATA_DIR=${MAIN_DIR}/data
RAW_READS_DIR=${DATA_DIR}/fastq
PRESORTED_READS_DIR=${DATA_DIR}/bam

RESULTS_DIR=${MAIN_DIR}/results
SSS_DIR=${RESULTS_DIR}/sss
SSS_SORTED_DIR=${SSS_DIR}/sorted_reads
SAMPLES_FILE=${RESULTS_DIR}/samples.tsv
LOG_FILE=${RESULTS_DIR}/log.txt

ENSEMBL_VERSION=86
MOUSE_GENOME_DIR=/srv/data/genome/mouse/ensembl-${ENSEMBL_VERSION}
MOUSE_STAR_INDEX=${MOUSE_GENOME_DIR}/STAR_indices/primary_assembly
RAT_GENOME_DIR=/srv/data/genome/rat/ensembl-${ENSEMBL_VERSION}
RAT_STAR_INDEX=${RAT_GENOME_DIR}/STAR_indices/toplevel
HUMAN_GENOME_DIR=/srv/data/genome/human/ensembl-${ENSEMBL_VERSION}
HUMAN_STAR_INDEX=${HUMAN_GENOME_DIR}/STAR_indices/primary_assembly

NUM_THREADS_PER_SAMPLE=1
NUM_TOTAL_THREADS=1
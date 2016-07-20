#!/bin/bash

function listFiles {
    DELIMITER=$1
    shift
    FILES=$@

    OUTPUT=`ls -1 $FILES | tr '\n' "$DELIMITER"`
    echo ${OUTPUT%$DELIMITER}
}

function create_samples_file_entry {
    SAMPLE=$1
    READS_FILES_1=$2
    READS_FILES_2=$3
    SAMPLES_FILE=$4

    echo "${SAMPLE}_reads ${READS_FILES_1} ${READS_FILES_2}" >> ${SAMPLES_FILE}
}

function merge_filtered_reads {
    NUM_THREADS=$1
    FILTERED_READS_DIR=$2
    SPECIES=$3
    SAMPLE=$4

    sambamba merge -t ${NUM_THREADS} ${FILTERED_READS_DIR}/${SAMPLE}_${SPECIES}.bam ${FILTERED_READS_DIR}/${SAMPLE}_reads_${SPECIES}*.bam
    rm ${FILTERED_READS_DIR}/${SAMPLE}_reads_${SPECIES}*filtered*
}

function count_reads_for_features {
    NUM_THREADS=$1
    FEATURES_GTF=$2
    BAM_FILE=$3
    COUNTS_OUTPUT_FILE=$4

    counts_tmp=.counts_tmp

    featureCounts -T ${NUM_THREADS} -p -a ${FEATURES_GTF} -o ${counts_tmp} ${BAM_FILE}
    tail -n +3 ${counts_tmp} | cut -f 1,7 > ${COUNTS_OUTPUT_FILE}

    rm ${counts_tmp}
    mv ${counts_tmp}.summary ${COUNTS_OUTPUT_FILE}.summary
}

function count_reads_for_features_stranded {
    NUM_THREADS=$1
    FEATURES_GTF=$2
    BAM_FILE=$3
    COUNTS_OUTPUT_FILE=$4

    counts_tmp=.counts_tmp

    featureCounts -T ${NUM_THREADS} -p -a ${FEATURES_GTF} -o ${counts_tmp} -s 2 ${BAM_FILE}
    tail -n +3 ${counts_tmp} | cut -f 1,7 > ${COUNTS_OUTPUT_FILE}

    rm ${counts_tmp}
    mv ${counts_tmp}.summary ${COUNTS_OUTPUT_FILE}.summary
}

function count_reads_in_gzipped_fastq_files {
    FASTQ_FILES=$@

    num_lines=$(zcat ${FASTQ_FILES} | wc -l)
    echo $[${num_lines}/4]
}

function sum_reads_mapped_to_features {
    COUNTS_FILE=$1

    echo $(awk '{s+=$2} END {print s}' ${COUNTS_FILE})
}

function tidy_differential_expression_results_file {
    RESULTS_FILE=$1

    sed -i "s/-Inf/'-Inf/g;s/,NA,/,,/g;s/,NA,/,,/g;s/,NA$/,/g" ${RESULTS_FILE}
}

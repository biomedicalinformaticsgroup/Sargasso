#!/bin/bash

set -o nounset
set -o errexit
set -o xtrace

source definitions.sh
source common_functions.sh

function run_sss_filter_conservative {
    for sample in ${MOUSE_ACTIVITY_SAMPLES}; do
        sample_reads_1=$(listFiles , ${MOUSE_ACTIVITY_RNA_SEQ_DIR}/${sample}/*_1.sanfastq.gz)
        sample_reads_2=$(listFiles , ${MOUSE_ACTIVITY_RNA_SEQ_DIR}/${sample}/*_2.sanfastq.gz)

        create_samples_file_entry ${sample} ${sample_reads_1} ${sample_reads_2} ${MOUSE_ACTIVITY_SAMPLES_FILE}
    done

    species_separator --conservative --reads-base="/" --s1-index=${MOUSE_STAR_INDEX} --s2-index=${RAT_STAR_INDEX} -t ${NUM_THREADS} mouse rat ${MOUSE_ACTIVITY_SAMPLES_FILE} ${MOUSE_ACTIVITY_SSS_DIR_CONSERVATIVE}

    (cd ${MOUSE_ACTIVITY_SSS_DIR_CONSERVATIVE}; make)
}

function run_sss_filter_best {
    species_separator --best --reads-base="/" --s1-index=${MOUSE_STAR_INDEX} --s2-index=${RAT_STAR_INDEX} -t ${NUM_THREADS} mouse rat ${MOUSE_ACTIVITY_SAMPLES_FILE} ${MOUSE_ACTIVITY_SSS_DIR_BEST}

    (cd ${MOUSE_ACTIVITY_SSS_DIR_BEST}; make)
}

function run_normal_mapping {
    species_separator --best --reads-base="/" --s1-index=${MOUSE_STAR_INDEX} --s2-index=${RAT_STAR_INDEX} -t ${NUM_THREADS} mouse rat ${MOUSE_ACTIVITY_SAMPLES_FILE} ${MOUSE_ACTIVITY_NORMAL_MAPPING_DIR}

    (cd ${MOUSE_ACTIVITY_NORMAL_MAPPING_DIR};
        make mapped_reads
    )
}

function count_reads {
    mkdir -p ${MOUSE_ACTIVITY_COUNTS_DIR}
    for sample in ${MOUSE_ACTIVITY_SAMPLES};
    do
        count_reads_for_features_stranded ${NUM_THREADS} ${MOUSE_GTF} ${MOUSE_ACTIVITY_SSS_FILTERED_READS_DIR_CONSERVATIVE}/${sample}_reads_mouse_filtered.bam ${MOUSE_ACTIVITY_COUNTS_DIR}/${sample}.conservative.count
        count_reads_for_features_stranded ${NUM_THREADS} ${MOUSE_GTF} ${MOUSE_ACTIVITY_SSS_FILTERED_READS_DIR_BEST}/${sample}_reads_mouse_filtered.bam ${MOUSE_ACTIVITY_COUNTS_DIR}/${sample}.best.count
        count_reads_for_features_stranded ${NUM_THREADS} ${MOUSE_GTF} ${MOUSE_ACTIVITY_NORMAL_MAPPING_MAPPED_READS_DIR}/${sample}_reads.mouse.bam ${MOUSE_ACTIVITY_COUNTS_DIR}/${sample}.normal.count
    done
}

function write_mapping_stats {
    stats_file=${MOUSE_ACTIVITY_RESULTS_DIR}/mapping_stats.csv

    echo "sample,read_pairs,normal,conservative,best" > ${stats_file}

    for sample in ${MOUSE_ACTIVITY_SAMPLES}; do
        echo -n "${sample}," >> ${stats_file}
        echo -n "$(count_reads_in_gzipped_fastq_files ${MOUSE_ACTIVITY_RNA_SEQ_DIR}/${sample}/*_1.sanfastq.gz)," >> ${stats_file}
        echo -n "$(sum_reads_mapped_to_features ${MOUSE_ACTIVITY_COUNTS_DIR}/${sample}.normal.count)," >> ${stats_file}
        echo "$(sum_reads_mapped_to_features ${MOUSE_ACTIVITY_COUNTS_DIR}/${sample}.conservative.count)" >> ${stats_file}
        echo "$(sum_reads_mapped_to_features ${MOUSE_ACTIVITY_COUNTS_DIR}/${sample}.best.count)" >> ${stats_file}
    done
}

function perform_differential_expression {
    diff_expr_dir=${MOUSE_ACTIVITY_RESULTS_DIR}/differential_expression
    mkdir -p ${diff_expr_dir}

    Rscript mouse_activity_differential_expression.R

    tidy_differential_expression_results_file ${diff_expr_dir}/diff_expr.csv
}

MOUSE_ACTIVITY_RESULTS_DIR=${RESULTS_DIR}/mouse_activity.test
MOUSE_ACTIVITY_RNA_SEQ_DIR=${RNASEQ_DIR}/neuron_astrocyte_activity_mouse
MOUSE_ACTIVITY_SSS_DIR_CONSERVATIVE=${MOUSE_ACTIVITY_RESULTS_DIR}/sss.conservative
MOUSE_ACTIVITY_SSS_DIR_BEST=${MOUSE_ACTIVITY_RESULTS_DIR}/sss.best
MOUSE_ACTIVITY_SSS_FILTERED_READS_DIR_CONSERVATIVE=${MOUSE_ACTIVITY_SSS_DIR_CONSERVATIVE}/filtered_reads
MOUSE_ACTIVITY_SSS_FILTERED_READS_DIR_BEST=${MOUSE_ACTIVITY_SSS_DIR_BEST}/filtered_reads
MOUSE_ACTIVITY_NORMAL_MAPPING_DIR=${MOUSE_ACTIVITY_RESULTS_DIR}/normal_mapping
MOUSE_ACTIVITY_NORMAL_MAPPING_MAPPED_READS_DIR=${MOUSE_ACTIVITY_NORMAL_MAPPING_DIR}/mapped_reads
MOUSE_ACTIVITY_COUNTS_DIR=${MOUSE_ACTIVITY_RESULTS_DIR}/read_counts

MOUSE_ACTIVITY_SAMPLES_FILE=${MOUSE_ACTIVITY_RESULTS_DIR}/samples.tsv
MOUSE_ACTIVITY_SAMPLES="AN1 AN2 BN1 BN2 CN1 CN2"

mkdir -p ${MOUSE_ACTIVITY_RESULTS_DIR}

run_sss_filter_conservative
run_sss_filter_best
run_normal_mapping
count_reads
write_mapping_stats
#perform_differential_expression

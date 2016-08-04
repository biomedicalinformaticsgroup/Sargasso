#!/bin/bash

set -o nounset
set -o errexit
set -o xtrace

source definitions.sh
source common_functions.sh

function get_protein_coding_exons {
    GTF_FILE=$1

    awk '$3 == "exon"' ${GTF_FILE} | grep  "gene_biotype \"protein_coding\""
}

function create_transcript_sequences {
    EXONS_GTF=$1
    CHR_FASTA=$2
    TRANSCRIPTS_PREFIX=$3

    rsem-prepare-reference --gtf ${EXONS_GTF} ${CHR_FASTA} ${TRANSCRIPTS_PREFIX}
    #rsem-prepare-reference --bowtie --gtf ${EXONS_GTF} ${CHR_FASTA} ${TRANSCRIPTS_PREFIX}
}

function create_reads_from_transcripts {
    READ_LENGTH=$1
    INSERT_SIZE=$2
    EXONS_GTF=$3
    TRANSCRIPTS_FASTA=$4
    READS_FILE_1=$5
    READS_FILE_2=$6

    READS_TMP=.tmp.fa

    transcripts_to_reads create --read-length ${READ_LENGTH} --insert-size ${INSERT_SIZE} --paired-end ${EXONS_GTF} ${TRANSCRIPTS_FASTA} > ${READS_TMP}

    grep -A 1 --no-group-separator "_1" ${READS_TMP} > ${READS_FILE_1}
    grep -A 1 --no-group-separator "_2" ${READS_TMP} > ${READS_FILE_2}

    gzip ${READS_FILE_1}
    gzip ${READS_FILE_2}

    rm ${READS_TMP}
}

function create_samples_file_entry {
    SAMPLE=$1
    READS_FILES_1=$2
    READS_FILES_2=$3
    SAMPLES_FILE=$4

    echo "${SAMPLE}_reads ${READS_FILES_1} ${READS_FILES_2}" >> ${SAMPLES_FILE}
}

function count_theoretical_reads_per_gene {
    READ_LENGTH=$1
    INSERT_SIZE=$2
    EXONS_GTF=$3
    TRANSCRIPTS_FASTA=$4

    transcripts_to_reads count --read-length ${READ_LENGTH} --insert-size ${INSERT_SIZE} --paired-end ${EXONS_GTF} ${TRANSCRIPTS_FASTA}
}

function count_mapped_reads {
    READ_COUNTS_FILE=$1
    BAM_FILE=$2
    MAPPING_TYPE=$3

    transcripts_to_reads count_mapped ${READ_COUNTS_FILE} ${BAM_FILE} |
        sed "s/mapped/${MAPPING_TYPE}_mapped/g"
}

function collate_mapped_read_counts {
    CORRECT_MAPPED=$1
    INCORRECT_MAPPED=$2

    join -t, ${CORRECT_MAPPED} ${INCORRECT_MAPPED} |
        cut -d, -f 1,2,3,4,6,7
}

function collate_assigned_read_counts {
    CONSERVATIVE_ASSIGNED=$1
    BEST_ASSIGNED=$2
    UNFILTERED_ASSIGNED=$3

    join ${CONSERVATIVE_ASSIGNED} ${BEST_ASSIGNED} | join - ${UNFILTERED_ASSIGNED}
}

#function collate_theoretical_data {
    #READ_COUNTS=$1
    #SPECIES=$2
    #READ_LENGTH=$3
    #OUTPUT_FILE=$4

    #Rscript collate_theoretical_data.R ${READ_COUNTS} ${SPECIES} ${READ_LENGTH} ${OUTPUT_FILE}
#}

function get_protein_coding_exons_for_both_species {
    get_protein_coding_exons ${MOUSE_GTF} > ${MOUSE_GTF_EXONS}
    get_protein_coding_exons ${RAT_GTF} > ${RAT_GTF_EXONS}
}

function create_transcript_sequences_for_both_species {
    create_transcript_sequences ${MOUSE_GTF_EXONS} ${MOUSE_CHR_FASTA} ${TRANSCRIPT_SEQUENCES_DIR}/${MOUSE_PREFIX}
    create_transcript_sequences ${RAT_GTF_EXONS} ${RAT_CHR_FASTA} ${TRANSCRIPT_SEQUENCES_DIR}/${RAT_PREFIX}
}

function create_reads_from_transcripts_for_both_species {
    create_reads_from_transcripts ${READ_LENGTH} ${INSERT_SIZE} ${MOUSE_GTF_EXONS} ${MOUSE_TRANSCRIPTS_FASTA} ${MOUSE_READS_1} ${MOUSE_READS_2}
    create_reads_from_transcripts ${READ_LENGTH} ${INSERT_SIZE} ${RAT_GTF_EXONS} ${RAT_TRANSCRIPTS_FASTA} ${RAT_READS_1} ${RAT_READS_2}
}

function separate_species_reads {
    samples_file=${THEORETICAL_RESULTS_DIR}/samples.tsv

    create_samples_file_entry mouse_reads ${MOUSE_READS_1}.gz ${MOUSE_READS_2}.gz ${samples_file}
    create_samples_file_entry rat_reads ${RAT_READS_1}.gz ${RAT_READS_2}.gz ${samples_file}

    species_separator --reads-base="/" --s1-index=${MOUSE_STAR_INDEX} --s2-index=${RAT_STAR_INDEX} -t ${NUM_THREADS} --conservative mouse rat ${samples_file} ${THEORETICAL_SSS_DIR_CONSERVATIVE}

    (cd ${THEORETICAL_SSS_DIR_CONSERVATIVE}; make)

    species_separator --reads-base="/" --s1-index=${MOUSE_STAR_INDEX} --s2-index=${RAT_STAR_INDEX} -t ${NUM_THREADS} --best mouse rat ${samples_file} ${THEORETICAL_SSS_DIR_BEST}

    (cd ${THEORETICAL_SSS_DIR_BEST}; make)
}

function count_theoretical_reads_per_gene_for_both_species {
    count_theoretical_reads_per_gene ${READ_LENGTH} ${INSERT_SIZE} ${MOUSE_GTF_EXONS} ${MOUSE_TRANSCRIPTS_FASTA} > ${MOUSE_READ_COUNTS}
    count_theoretical_reads_per_gene ${READ_LENGTH} ${INSERT_SIZE} ${RAT_GTF_EXONS} ${RAT_TRANSCRIPTS_FASTA} > ${RAT_READ_COUNTS}
}

function count_mapped_reads_for_both_species {
    count_mapped_reads ${MOUSE_READ_COUNTS} ${THEORETICAL_SSS_FILTERED_READS_DIR_CONSERVATIVE}/mouse_reads_reads_mouse_filtered.bam correct > ${MOUSE_CORRECT_READS_PER_GENE_CONSERVATIVE}
    count_mapped_reads ${RAT_READ_COUNTS} ${THEORETICAL_SSS_FILTERED_READS_DIR_CONSERVATIVE}/rat_reads_reads_rat_filtered.bam correct > ${RAT_CORRECT_READS_PER_GENE_CONSERVATIVE}
    count_mapped_reads ${MOUSE_READ_COUNTS} ${THEORETICAL_SSS_FILTERED_READS_DIR_CONSERVATIVE}/mouse_reads_reads_rat_filtered.bam incorrect > ${MOUSE_INCORRECT_READS_PER_GENE_CONSERVATIVE}
    count_mapped_reads ${RAT_READ_COUNTS} ${THEORETICAL_SSS_FILTERED_READS_DIR_CONSERVATIVE}/rat_reads_reads_mouse_filtered.bam incorrect > ${RAT_INCORRECT_READS_PER_GENE_CONSERVATIVE}

    count_mapped_reads ${MOUSE_READ_COUNTS} ${THEORETICAL_SSS_FILTERED_READS_DIR_BEST}/mouse_reads_reads_mouse_filtered.bam correct > ${MOUSE_CORRECT_READS_PER_GENE_BEST}
    count_mapped_reads ${RAT_READ_COUNTS} ${THEORETICAL_SSS_FILTERED_READS_DIR_BEST}/rat_reads_reads_rat_filtered.bam correct > ${RAT_CORRECT_READS_PER_GENE_BEST}
    count_mapped_reads ${MOUSE_READ_COUNTS} ${THEORETICAL_SSS_FILTERED_READS_DIR_BEST}/mouse_reads_reads_rat_filtered.bam incorrect > ${MOUSE_INCORRECT_READS_PER_GENE_BEST}
    count_mapped_reads ${RAT_READ_COUNTS} ${THEORETICAL_SSS_FILTERED_READS_DIR_BEST}/rat_reads_reads_mouse_filtered.bam incorrect > ${RAT_INCORRECT_READS_PER_GENE_BEST}
}

function collate_mapped_read_counts_for_both_species {
    collate_mapped_read_counts ${MOUSE_CORRECT_READS_PER_GENE_CONSERVATIVE} ${MOUSE_INCORRECT_READS_PER_GENE_CONSERVATIVE} > ${MOUSE_READS_PER_GENE_CONSERVATIVE}
    collate_mapped_read_counts ${RAT_CORRECT_READS_PER_GENE_CONSERVATIVE} ${RAT_INCORRECT_READS_PER_GENE_CONSERVATIVE} > ${RAT_READS_PER_GENE_CONSERVATIVE}

    collate_mapped_read_counts ${MOUSE_CORRECT_READS_PER_GENE_BEST} ${MOUSE_INCORRECT_READS_PER_GENE_BEST} > ${MOUSE_READS_PER_GENE_BEST}
    collate_mapped_read_counts ${RAT_CORRECT_READS_PER_GENE_BEST} ${RAT_INCORRECT_READS_PER_GENE_BEST} > ${RAT_READS_PER_GENE_BEST}
}

function count_assigned_reads_for_both_species {
    count_reads_for_features ${NUM_THREADS} ${MOUSE_GTF_EXONS} ${THEORETICAL_SSS_FILTERED_READS_DIR_CONSERVATIVE}/mouse_reads_reads_mouse_filtered.bam ${MOUSE_ASSIGNED_READS_PER_GENE_CONSERVATIVE}
    count_reads_for_features ${NUM_THREADS} ${MOUSE_GTF_EXONS} ${THEORETICAL_SSS_FILTERED_READS_DIR_BEST}/mouse_reads_reads_mouse_filtered.bam ${MOUSE_ASSIGNED_READS_PER_GENE_BEST}
    count_reads_for_features ${NUM_THREADS} ${MOUSE_GTF_EXONS} ${THEORETICAL_SSS_MAPPED_READS_DIR_BEST}/mouse_reads_reads.mouse.bam ${MOUSE_ASSIGNED_READS_PER_GENE_UNFILTERED}

    count_reads_for_features ${NUM_THREADS} ${RAT_GTF_EXONS} ${THEORETICAL_SSS_FILTERED_READS_DIR_CONSERVATIVE}/rat_reads_reads_rat_filtered.bam ${RAT_ASSIGNED_READS_PER_GENE_CONSERVATIVE}
    count_reads_for_features ${NUM_THREADS} ${RAT_GTF_EXONS} ${THEORETICAL_SSS_FILTERED_READS_DIR_BEST}/rat_reads_reads_rat_filtered.bam ${RAT_ASSIGNED_READS_PER_GENE_BEST}
    count_reads_for_features ${NUM_THREADS} ${RAT_GTF_EXONS} ${THEORETICAL_SSS_MAPPED_READS_DIR_BEST}/rat_reads_reads.rat.bam ${RAT_ASSIGNED_READS_PER_GENE_UNFILTERED}
}

function collate_assigned_read_counts_for_both_species {
    collate_assigned_read_counts ${MOUSE_ASSIGNED_READS_PER_GENE_CONSERVATIVE} ${MOUSE_ASSIGNED_READS_PER_GENE_BEST} ${MOUSE_ASSIGNED_READS_PER_GENE_UNFILTERED} > ${MOUSE_ASSIGNED_READS_PER_GENE}
    collate_assigned_read_counts ${RAT_ASSIGNED_READS_PER_GENE_CONSERVATIVE} ${RAT_ASSIGNED_READS_PER_GENE_BEST} ${RAT_ASSIGNED_READS_PER_GENE_UNFILTERED} > ${RAT_ASSIGNED_READS_PER_GENE}
}

#function collate_theoretical_data_for_both_species {
    #collate_theoretical_data ${MOUSE_READS_PER_GENE_CONSERVATIVE} mouse ${READ_LENGTH} ${THEORETICAL_RESULTS_DIR}/mouse_feasibility.conservative.csv
    #collate_theoretical_data ${RAT_READS_PER_GENE_CONSERVATIVE} rat ${READ_LENGTH} ${THEORETICAL_RESULTS_DIR}/rat_feasibility.conservative.csv

    #collate_theoretical_data ${MOUSE_READS_PER_GENE_BEST} mouse ${READ_LENGTH} ${THEORETICAL_RESULTS_DIR}/mouse_feasibility.best.csv
    #collate_theoretical_data ${RAT_READS_PER_GENE_BEST} rat ${READ_LENGTH} ${THEORETICAL_RESULTS_DIR}/rat_feasibility.best.csv
#}

THEORETICAL_RESULTS_DIR=${RESULTS_DIR}/theoretical
TRANSCRIPT_SEQUENCES_DIR=${THEORETICAL_RESULTS_DIR}/transcript_sequences
THEORETICAL_SSS_DIR_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/sss.conservative
THEORETICAL_SSS_FILTERED_READS_DIR_CONSERVATIVE=${THEORETICAL_SSS_DIR_CONSERVATIVE}/filtered_reads
THEORETICAL_SSS_DIR_BEST=${THEORETICAL_RESULTS_DIR}/sss.best
THEORETICAL_SSS_FILTERED_READS_DIR_BEST=${THEORETICAL_SSS_DIR_BEST}/filtered_reads
THEORETICAL_SSS_MAPPED_READS_DIR_BEST=${THEORETICAL_SSS_DIR_BEST}/mapped_reads

READ_LENGTH=50
INSERT_SIZE=150

MOUSE_PREFIX=mm_ensembl${ENSEMBL_VERSION}
RAT_PREFIX=rn_ensembl${ENSEMBL_VERSION}

MOUSE_GTF_EXONS=${TRANSCRIPT_SEQUENCES_DIR}/${MOUSE_PREFIX}.exons.gtf
RAT_GTF_EXONS=${TRANSCRIPT_SEQUENCES_DIR}/${RAT_PREFIX}.exons.gtf

MOUSE_READS_1=${TRANSCRIPT_SEQUENCES_DIR}/${MOUSE_PREFIX}.reads_${READ_LENGTH}_1.fa
MOUSE_READS_2=${TRANSCRIPT_SEQUENCES_DIR}/${MOUSE_PREFIX}.reads_${READ_LENGTH}_2.fa

RAT_READS_1=${TRANSCRIPT_SEQUENCES_DIR}/${RAT_PREFIX}.reads_${READ_LENGTH}_1.fa
RAT_READS_2=${TRANSCRIPT_SEQUENCES_DIR}/${RAT_PREFIX}.reads_${READ_LENGTH}_2.fa

MOUSE_TRANSCRIPTS_FASTA=${TRANSCRIPT_SEQUENCES_DIR}/${MOUSE_PREFIX}.transcripts.fa
RAT_TRANSCRIPTS_FASTA=${TRANSCRIPT_SEQUENCES_DIR}/${RAT_PREFIX}.transcripts.fa

MOUSE_READ_COUNTS=${TRANSCRIPT_SEQUENCES_DIR}/${MOUSE_PREFIX}.read_counts_${READ_LENGTH}.csv
RAT_READ_COUNTS=${TRANSCRIPT_SEQUENCES_DIR}/${RAT_PREFIX}.read_counts_${READ_LENGTH}.csv

MOUSE_CORRECT_READS_PER_GENE_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/unambiguous_mouse_reads_${READ_LENGTH}.conservative.csv
RAT_CORRECT_READS_PER_GENE_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/unambiguous_rat_reads_${READ_LENGTH}.conservative.csv
MOUSE_CORRECT_READS_PER_GENE_BEST=${THEORETICAL_RESULTS_DIR}/unambiguous_mouse_reads_${READ_LENGTH}.best.csv
RAT_CORRECT_READS_PER_GENE_BEST=${THEORETICAL_RESULTS_DIR}/unambiguous_rat_reads_${READ_LENGTH}.best.csv

MOUSE_INCORRECT_READS_PER_GENE_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/incorrect_mouse_reads_${READ_LENGTH}.conservative.csv
RAT_INCORRECT_READS_PER_GENE_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/incorrect_rat_reads_${READ_LENGTH}.conservative.csv
MOUSE_INCORRECT_READS_PER_GENE_BEST=${THEORETICAL_RESULTS_DIR}/incorrect_mouse_reads_${READ_LENGTH}.best.csv
RAT_INCORRECT_READS_PER_GENE_BEST=${THEORETICAL_RESULTS_DIR}/incorrect_rat_reads_${READ_LENGTH}.best.csv

MOUSE_ASSIGNED_READS_PER_GENE_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/assigned_mouse_reads_${READ_LENGTH}.conservative.csv
MOUSE_ASSIGNED_READS_PER_GENE_BEST=${THEORETICAL_RESULTS_DIR}/assigned_mouse_reads_${READ_LENGTH}.best.csv
MOUSE_ASSIGNED_READS_PER_GENE_UNFILTERED=${THEORETICAL_RESULTS_DIR}/assigned_mouse_reads_${READ_LENGTH}.unfiltered.csv

RAT_ASSIGNED_READS_PER_GENE_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/assigned_rat_reads_${READ_LENGTH}.conservative.csv
RAT_ASSIGNED_READS_PER_GENE_BEST=${THEORETICAL_RESULTS_DIR}/assigned_rat_reads_${READ_LENGTH}.best.csv
RAT_ASSIGNED_READS_PER_GENE_UNFILTERED=${THEORETICAL_RESULTS_DIR}/assigned_rat_reads_${READ_LENGTH}.unfiltered.csv

MOUSE_READS_PER_GENE_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/mouse_reads_per_gene_${READ_LENGTH}.conservative.csv 
RAT_READS_PER_GENE_CONSERVATIVE=${THEORETICAL_RESULTS_DIR}/rat_reads_per_gene_${READ_LENGTH}.conservative.csv 
MOUSE_READS_PER_GENE_BEST=${THEORETICAL_RESULTS_DIR}/mouse_reads_per_gene_${READ_LENGTH}.best.csv 
RAT_READS_PER_GENE_BEST=${THEORETICAL_RESULTS_DIR}/rat_reads_per_gene_${READ_LENGTH}.best.csv 

MOUSE_ASSIGNED_READS_PER_GENE=${THEORETICAL_RESULTS_DIR}/assigned_mouse_reads_per_gene_${READ_LENGTH}.csv
RAT_ASSIGNED_READS_PER_GENE=${THEORETICAL_RESULTS_DIR}/assigned_rat_reads_per_gene_${READ_LENGTH}.csv

mkdir -p ${TRANSCRIPT_SEQUENCES_DIR}

get_protein_coding_exons_for_both_species
create_transcript_sequences_for_both_species
create_reads_from_transcripts_for_both_species
separate_species_reads
count_theoretical_reads_per_gene_for_both_species
count_mapped_reads_for_both_species
collate_mapped_read_counts_for_both_species
count_assigned_reads_for_both_species
collate_assigned_read_counts_for_both_species

#collate_theoretical_data_for_both_species

#!/bin/bash

NUM_THREADS=16

# EDIT THIS
MAIN_DIR=/home/odando/projects/new_species_separator/paper_analysis/

DATA_DIR=${MAIN_DIR}/data
# TODO: setup links to mouse and rat Ensembl directories

RNASEQ_DIR=${DATA_DIR}/rnaseq
RESULTS_DIR=${MAIN_DIR}/results

ENSEMBL_VERSION=84
MOUSE_GENOME_DIR=${DATA_DIR}/mouse_ensembl-${ENSEMBL_VERSION}
RAT_GENOME_DIR=${DATA_DIR}/rat_ensembl-${ENSEMBL_VERSION}

MOUSE_GTF=${MOUSE_GENOME_DIR}/Mus_musculus.GRCm38.${ENSEMBL_VERSION}.gtf
MOUSE_CHR_FASTA=${MOUSE_GENOME_DIR}/primary_assembly
MOUSE_STAR_INDEX=${MOUSE_GENOME_DIR}/STAR_indices/primary_assembly

RAT_GTF=${RAT_GENOME_DIR}/Rattus_norvegicus.Rnor_6.0.${ENSEMBL_VERSION}.gtf
RAT_CHR_FASTA=${RAT_GENOME_DIR}/toplevel
RAT_STAR_INDEX=${RAT_GENOME_DIR}/STAR_indices/toplevel


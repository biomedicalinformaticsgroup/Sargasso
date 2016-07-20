#!/bin/bash

source definitions.sh

#./setup.sh

#get_gene_lengths ${MOUSE_GTF} > ${RESULTS_DIR}/gene_lengths.csv
#get_gene_lengths ${RAT_GTF} > ${RESULTS_DIR}/rat_gene_lengths.csv

./create_theoretical_results.sh

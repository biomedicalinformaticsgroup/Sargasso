#!/usr/bin/env bash

# to test chipseq PE
species_separator chipseq --run-separation \
    --reads-base-dir=/home/xinhe/Projects/Sargasso/pipeline_test/data/fastq/ \
    --num-threads 4 --best \
    --sambamba-sort-tmp-dir /home/xinhe/tmp /home/xinhe/Projects/Sargasso/pipeline_test/data/chipseq.tsv \
    /home/xinhe/Projects/Sargasso/results/chipseqtest \
    mouse /srv/data/genome/mouse/ensembl-93/bowtie2_indexes/primary_assembly \
    human /srv/data/genome/human/ensembl-93/bowtie2_indexes/primary_assembly \
    rat /srv/data/genome/rat/ensembl-93/bowtie2_indexes/toplevel


# to test chipseq PE + index build
species_separator chipseq --run-separation \
    --reads-base-dir=/home/xinhe/Projects/Sargasso/pipeline_test/data/fastq/ \
    --num-threads 4 --best \
    --sambamba-sort-tmp-dir /home/xinhe/tmp /home/xinhe/Projects/Sargasso/pipeline_test/data/chipseq.tsv \
    /home/xinhe/Projects/Sargasso/results/chipseqtest \
    mouse /srv/data/genome/mouse/ensembl-93/mouse_primary_assembly.fa \
    human /srv/data/genome/human/ensembl-93/human_primary_assembly.fa \
    rat /srv/data/genome/rat/ensembl-93/rat_toplevel.fa

# to test chipseq SE
species_separator chipseq --run-separation \
    --reads-base-dir=/home/xinhe/Projects/Sargasso/pipeline_test/data/fastq/ \
    --num-threads 4 --best \
    --sambamba-sort-tmp-dir /home/xinhe/tmp /home/xinhe/Projects/Sargasso/pipeline_test/data/chipseq_SE.tsv \
    /home/xinhe/Projects/Sargasso/results/chipseqtest \
    mouse /srv/data/genome/mouse/ensembl-93/mouse_primary_assembly.fa \
    human /srv/data/genome/human/ensembl-93/human_primary_assembly.fa \
    rat /srv/data/genome/rat/ensembl-93/rat_toplevel.fa
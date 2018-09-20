#!/usr/bin/env bash

#https://www.ncbi.nlm.nih.gov/gds
#
#with search filter:
#("neurology"[MeSH Terms] OR neuro[All Fields]) AND "Genome binding/occupancy profiling by high throughput sequen"[Filter]
#
#https://www.ncbi.nlm.nih.gov/gds?term=%28%22neurology%22%5BMeSH%20Terms%5D%20OR%20neuro%5BAll%20Fields%5D%29%20AND%20%22Genome%20binding/occupancy%20profiling%20by%20high%20throughput%20sequen%22%5BFilter%5D&cmd=DetailsSearch
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75942


#paired end example
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110032
#ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR668/SRR6685151/SRR6685151.sra
#fastq-dump --split-3 SRR6685151.sra &
#gzip chipseq_mouse_R1.fastq
#gzip chipseq_mouse_R2.fastq




#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR298/SRR2989997/SRR2989997.sra
#fastq-dump --split-3 SRR2989997.sra &
#head -400 SRR2989997.fastq > chipseq_mouse.fastq
#gzip chipseq_mouse.fastq
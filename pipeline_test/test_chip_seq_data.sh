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
head -40000 SRR6685151_1.fastq > chipseq_mouse_R1.fastq
head -40000 SRR6685151_2.fastq  > chipseq_mouse_R2.fastq
gzip chipseq_mouse_R1.fastq
gzip chipseq_mouse_R2.fastq







#single end example
#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR298/SRR2989997/SRR2989997.sra
#fastq-dump --split-3 SRR2989997.sra &
#head -400 SRR2989997.fastq > chipseq_mouse.fastq
#gzip chipseq_mouse.fastq


bowtie2-build --threads 48 /srv/data/genome/rat/ensembl-93/rat_toplevel.fa rat93
bowtie2-build --threads 48 /srv/data/genome/mouse/ensembl-93/mouse_primary_assembly.fa  mouse93
bowtie2-build --threads 48 /srv/data/genome/human/ensembl-93/human_primary_assembly.fa  human93

bowtie2 --no-unal --no-discordant --no-mixed -k 100 --threads 24  -x rat93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S rat_paired.sam  >rat_paired.log 2>&1 &
bowtie2 --no-unal --no-discordant --no-mixed -k 100 --threads 24  -x mouse93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S mouse_paired.sam >mouse_paired.log 2>&1 &

sambamba view -S mouse_paired.sam -f bam > mouse_paired.bam
sambamba view -S rat_paired.sam -f bam > rat_paired.bam


mv mouse_paired.bam chiseq_mouse_sample.mouse.bam
mv rat_paired.bam chiseq_mouse_sample.rat.bam
sort_reads "mouse rat" chiseq_mouse_sample 2 /home/xinhe/Projects/Sargasso/results/chipseq/mapped_reads /home/xinhe/Projects/Sargasso/results/chipseq/sorted_reads ~/tmp
filter_reads rnaseq "chiseq_mouse_sample" /home/xinhe/Projects/Sargasso/results/chipseq/sorted_reads /home/xinhe/Projects/Sargasso/results/chipseq/filtered_reads 2 1 2 999999 "" mouse rat

filter_control rnaseq /home/xinhe/Projects/Sargasso/results/chipseq/filtered_reads/Blocks /home/xinhe/Projects/Sargasso/results/chipseq/filtered_reads chiseq_mouse_sample 1 2 999999 mouse rat


/usr/bin/time bowtie2 -x mouse93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S mouse_paired_default.sam
/usr/bin/time bowtie2 -k 1 -x mouse93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S mouse_paired_1.sam
/usr/bin/time bowtie2 --no-unal --no-discordant --no-mixed -p 24 -k 100 -x mouse93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S mouse_paired_100.sam
/usr/bin/time bowtie2 -k 10000 -x mouse93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S mouse_paired_10000.sam


 samtools view -F 4 -S mouse_paired_100.sam | grep -v ^@ | cut -f 1 | uniq -c | sort -n -r -k 1 | cut -f 1 -d'S' | uniq -c



bowtie2 --no-unal --no-discordant --no-mixed -p 24 -k 100 -x mouse93 -1 chipseq_rat_R1.fastq.gz -2 chipseq_rat_R2.fastq.gz -S rat_paired.sam
bowtie2 --no-unal --no-discordant --no-mixed -p 24 -k 100 -x mouse93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S mouse_paired.sam
sambamba view -S mouse_paired.sam -f bam > chiseq_mouse_sample.bam
sambamba view -S rat_paired.sam -f bam > chiseq_mouse_sample.bam
cp chiseq_mouse_sample.mouse.bam chiseq_mouse_sample.rat.bam /home/xinhe/Projects/Sargasso/results/chipseq/mapped_reads
sort_reads "mouse rat" chiseq_mouse_sample 2 /home/xinhe/Projects/Sargasso/results/chipseq/mapped_reads /home/xinhe/Projects/Sargasso/results/chipseq/sorted_reads ~/tmp
filter_reads rnaseq "chiseq_mouse_sample" sorted_reads filtered_reads 1 3 2 999999 "" mouse rat
multiqc -d -f -m sargasso .



head -4000000 SRR6685151_1.fastq > chipseq_mouse_R1.fastq &
head -4000000 SRR6685151_2.fastq  > chipseq_mouse_R2.fastq &
gzip chipseq_mouse_R1.fastq
gzip chipseq_mouse_R2.fastq

bowtie2 --no-unal --no-discordant --no-mixed -k 100 --threads 24  -x rat93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S rat_paired.sam  >rat_paired.log 2>&1 &
bowtie2 --no-unal --no-discordant --no-mixed -k 100 --threads 24  -x mouse93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S mouse_paired.sam >mouse_paired.log 2>&1 &
bowtie2 --no-unal --no-discordant --no-mixed -k 100 --threads 24  -x human93 -1 chipseq_mouse_R1.fastq.gz -2 chipseq_mouse_R2.fastq.gz -S human_paired.sam >human_paired.log 2>&1 &

sambamba view -S mouse_paired.sam -f bam > /home/xinhe/Projects/Sargasso/results/chipseq/mapped_reads/chiseq_mouse_sample.mouse.bam &
sambamba view -S rat_paired.sam -f bam > /home/xinhe/Projects/Sargasso/results/chipseq/mapped_reads/chiseq_mouse_sample.rat.bam &
sambamba view -S human_paired.sam -f bam > /home/xinhe/Projects/Sargasso/results/chipseq/mapped_reads/chiseq_mouse_sample.human.bam &

sort_reads "mouse rat human" chiseq_mouse_sample 24 /home/xinhe/Projects/Sargasso/results/chipseq/mapped_reads /home/xinhe/Projects/Sargasso/results/chipseq/sorted_reads ~/tmp
filter_reads rnaseq "chiseq_mouse_sample" sorted_reads filtered_reads 1 1 2 999999 "" mouse rat human



filter_reads rnaseq "chiseq_mouse_sample" sorted_reads filtered_reads 1 1 2 999999 "" mouse rat human

multiqc -d -f -m sargasso -o filtered_reads_human_mouse_rat filtered_reads_human_mouse_rat
#!/usr/bin/env bash
#bash /home/xinhe/Projects/Sargasso/bin/sargasso_parameter_test.sh --output_dir ~/tmp/sargasso_test2 \
#--rnaseq_dir '/srv/data/ghardingham/snakemake_test' \
#--read1_identifier '01_1' \
#--read2_identifier '01_2' \
#--fastq_suffix 'fastq.gz' \
#--samples '1467_A 1467_V 1468_P' \
#--samples_origin 'mouse human rat' \
#--mismatch_setting '0 2 4 ' \
#--minmatch_setting '0 2 5' \
#--mutlimap_setting '1' \
#--mapper_executable STAR2.7.0f \
#--num_total_threads 16 \
#--species_para 'human /srv/data/genome/human/ensembl-99/STAR_indices/primary_assembly' \
#--species_para 'mouse /srv/data/genome/mouse/ensembl-99/STAR_indices/primary_assembly' \
#--species_para 'rat /srv/data/genome/rat/ensembl-99/STAR_indices/toplevel'

function usage {
  cat <<EOT

bash ./sargasso_parameter_test.sh --output_dir ~/tmp/sargasso_test \
--rnaseq_dir '/srv/data/ghardingham/snakemake_test' \
--read1_identifier '01_1' \
--read2_identifier '01_2' \
--fastq_suffix 'fastq.gz' \
--paired_end_read 'no' \
--samples '1467_A 1467_V 1468_P' \
--samples_origin 'mouse mouse rat' \
--mismatch_setting '0 2 4' \
--minmatch_setting '0 2 5' \
--mutlimap_setting '1' \
--mapper_executable STAR2.7.0f \
--num_total_threads 16 \
--plot_format pdf \
--species_para 'human /srv/data/genome/human/ensembl-99/STAR_indices/primary_assembly' \
--species_para 'mouse /srv/data/genome/mouse/ensembl-99/STAR_indices/primary_assembly' \
--species_para 'rat /srv/data/genome/rat/ensembl-99/STAR_indices/toplevel'

Usage:
  $(basename $0)
    [--help]
    --output_dir /tmp/sargasso_test
    --rnaseq_dir '/srv/data/ghardingham/snakemake_test'
    --samples '1467_A 1467_V 1468_P'
    --samples_origin 'mouse mouse rat'
    --species_para 'human /srv/data/genome/human/ensembl-99/STAR_indices/primary_assembly'
    --species_para 'mouse /srv/data/genome/mouse/ensembl-99/STAR_indices/primary_assembly'
    --species_para 'rat /srv/data/genome/rat/ensembl-99/STAR_indices/toplevel'
    [--read1_identifier '01_1']
    [--read2_identifier '01_2']
    [--fastq_suffix 'fastq.gz']
    [--paired_end_read]
    [--skip_init_run]
    [--mismatch_setting '0 2 4']
    [--minmatch_setting '0 2 5']
    [--mutlimap_setting '1']
    [--mapper_executable=STAR]
    [--num_total_threads=4]
    [--plot_format pdf]

--help
        see this help message

--output_dir
        all output file will be stored in this folder.

--samples
        quoted space separated sample names. Example "sample1 sample2 sample3"

--samples_origin
        quoted space separated sample original species. Example "mouse mouse rat"

--species_para
        quoted space separated species/genome index pair. Example: 'human /srv/data/genome/human/ensembl-99/STAR_indices/primary_assembly'

        this parameter can be provided multiple times to test sargasso on multiple species
        Example:        --species_para 'human genome/human/ensembl-99/STAR_indices/primary_assembly'
                        --species_para 'mouse genome/mouse/ensembl-99/STAR_indices/primary_assembly'
                        --species_para 'rat genome/rat/ensembl-99/STAR_indices/toplevel'

--read1_identifier
        pattern to identify raw read file. [default 001_1]

--read2_identifier
        pattern to identify raw read file for the second pair of reads in a paird_end_read data set. [default 001_2]

--fastq_suffix
        [default fastq.gz]

--paired_end_read
        is the data pair_end_read

--skip_init_run
        with multiple parameter to be tested, the script only create one initial mapping to save time
        if provide, the script will skip the initial run and directly perform the parameter test

--mismatch_setting
--minmatch_setting
--mutlimap_setting
        quoted space separated numbers to be tested. Example: "0 2 4 6"

--mapper_executable
        [default STAR]

--num_total_threads
        [default 1]

--plot_format
        can be one of pdf, png
        [default pdf]


EOT
}

set -o nounset
#set -o errexit
#set -o xtrace


OPTS="$(getopt -o h -l help,output_dir:,rnaseq_dir:,samples:,samples_origin:,read1_identifier:,read2_identifier:,fastq_suffix:,paired_end_read,skip_init_run,mismatch_setting:,minmatch_setting:,mutlimap_setting:,mapper_executable:,num_total_threads:,plot_format:,species_para: --name "$(basename "$0")" -- "$@")"

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi
if [ "$#" -eq  "0"  ] ; then usage ; exit 1 ; fi



READ1_IDENTIFIER='01_1'
READ2_IDENTIFIER='01_2'
FASTQ_SUFFIX='fastq.gz'
PAIRED_END_READ='no'
SKIP_INIT_RUN='no'
MISMATCH_SETTING='0 2 4 6 8 10'
MINMATCH_SETTING='0 2 5 10 '
MUTLIMAP_SETTING='1'
MAPPER_EXECUTABLE=STAR2.7.0f
NUM_TOTAL_THREADS=1
PLOT_FORMAT=pdf
SPECIES_PARA=()

eval set -- "$OPTS"
while true; do
  case $1 in
    -h | --help ) usage;;
    --output_dir ) OUTPUT_DIR=$2; shift 2 ;;
    --rnaseq_dir ) RNASEQ_DIR=$2; shift 2 ;;
    --samples ) SAMPLES=$2; shift 2 ;;
    --samples_origin ) SAMPLES_ORIGIN=$2; shift 2 ;;
    --read1_identifier ) READ1_IDENTIFIER=$2; shift 2 ;;
    --read2_identifier ) READ2_IDENTIFIER=$2; shift 2 ;;
    --fastq_suffix ) FASTQ_SUFFIX=$2; shift 2 ;;
    --paired_end_read ) PAIRED_END_READ='yes'; shift ;;
    --skip_init_run ) SKIP_INIT_RUN='yes'; shift ;;
    --mismatch_setting ) MISMATCH_SETTING=$2; shift 2 ;;
    --minmatch_setting ) MINMATCH_SETTING=$2; shift 2 ;;
    --mutlimap_setting ) MUTLIMAP_SETTING=$2; shift 2 ;;
    --mapper_executable ) MAPPER_EXECUTABLE=$2; shift 2 ;;
    --num_total_threads ) NUM_TOTAL_THREADS=$2; shift 2 ;;
    --plot_format ) PLOT_FORMAT=$2; shift 2 ;;
    --species_para ) SPECIES_PARA+=($2); shift 2 ;;
    -- ) shift; break  ;;
    * ) usage ;;
  esac
done

#echo OUTPUT_DIR: $OUTPUT_DIR
#echo RNASEQ_DIR: $RNASEQ_DIR
#echo SAMPLES: $SAMPLES
#echo SAMPLES_ORIGIN: $SAMPLES_ORIGIN
#echo READ1_IDENTIFIER: $READ1_IDENTIFIER
#echo READ2_IDENTIFIER: $READ2_IDENTIFIER
#echo FASTQ_SUFFIX: $FASTQ_SUFFIX
#echo PAIRED_END_READ: $PAIRED_END_READ
#echo SKIP_INIT_RUN: SKIP_INIT_RUN
#echo MISMATCH_SETTING: $MISMATCH_SETTING
#echo MINMATCH_SETTING: $MINMATCH_SETTING
#echo MUTLIMAP_SETTING: $MUTLIMAP_SETTING
#echo MAPPER_EXECUTABLE: $MAPPER_EXECUTABLE
#echo NUM_TOTAL_THREADS: $NUM_TOTAL_THREADS
#echo SPECIES_PARA: ${SPECIES_PARA[@]}
#echo PLOT_FORMAT: ${PLOT_FORMAT}


## we check mandatory parameters
[[ -z ${OUTPUT_DIR} ]] && echo "Error: --output_dir is missing!" && exit 1
[[ -z ${RNASEQ_DIR} ]] && echo "Error: --rnaseq_dir is missing!" && exit 1
[[ -z ${SAMPLES} ]] && echo "Error: --samples is missing!" && exit 1
[[ -z ${SAMPLES_ORIGIN} ]] && echo "Error: --samples_origin is missing!" && exit 1
[[ -z ${SPECIES_PARA} ]] && echo "Error: --species_para is missing!" && exit 1

## we check --samples" and "--samples_origin" have the same number of elements
number_sample=`echo "${SAMPLES}" | awk -F' ' '{print NF}'`
number_sample_origin=`echo "${SAMPLES_ORIGIN}" | awk -F' ' '{print NF}'`
[[ ${number_sample} -ne ${number_sample} ]] && echo "Error: number of sample does not equal to number of sample origin." && exit 1

## we only support png and pdf for now
[[ ! "${PLOT_FORMAT}" =~ ^(pdf|png)$ ]] && echo "Error: --plot_format can only be one of pdf/png!" && exit 1


mkdir -p ${OUTPUT_DIR}

function listFilesNoNewLine {
    local DELIMITER=$1
    shift
    local FILES=$@

    local OUTPUT=$(ls -1 $FILES | tr '\n' "${DELIMITER}")
    echo -n ${OUTPUT%$DELIMITER}
}

function addSample2tsv {
    local SAMPLE_TSV=$1 && shift
    local BASE_DIR=$1 && shift
    local READ1_IDENTIFIER=$1 && shift
    local READ2_IDENTIFIER=$1 && shift
    local FASTQ_SUFFIX=$1 && shift
    local PAIRED_READ=$1 && shift
    if [ -z "$READ2_IDENTIFIER" ];then
        local PAIRED_READ=0
    fi
    echo -ne '' > ${SAMPLE_TSV}
    local SAMPLE=$@
    for sample in ${SAMPLE}; do
        echo -ne ${sample}" " >> ${SAMPLE_TSV}
        echo -n $(listFilesNoNewLine "," ${BASE_DIR}/${sample}/*${READ1_IDENTIFIER}.${FASTQ_SUFFIX}) >> ${SAMPLE_TSV}
        if [ "$PAIRED_READ" == "yes" ]; then
            echo -n " "  >> ${SAMPLE_TSV}
            echo -n $(listFilesNoNewLine "," ${BASE_DIR}/${sample}/*${READ2_IDENTIFIER}.${FASTQ_SUFFIX}) >> ${SAMPLE_TSV}
        fi
        echo "" >> ${SAMPLE_TSV}
    done
}




TMP_DIR=${OUTPUT_DIR}/tmp
init_dir=${OUTPUT_DIR}/init_run
mkdir -p  ${TMP_DIR}
SAMPLE_TSV=${OUTPUT_DIR}/sample.tsv
addSample2tsv ${SAMPLE_TSV} ${RNASEQ_DIR} ${READ1_IDENTIFIER} ${READ2_IDENTIFIER} ${FASTQ_SUFFIX} ${PAIRED_END_READ} ${SAMPLES}


if [  "${SKIP_INIT_RUN}" == "no"  ]; then
    echo "Running Sargasso initial mapping...."
    ## we run the first sargasso run to get the mapped reads
    species_separator rnaseq --mapper-executable ${MAPPER_EXECUTABLE}  --sambamba-sort-tmp-dir=${TMP_DIR} \
                    --mismatch-threshold 0 --minmatch-threshold 0 --multimap-threshold 1 --reject-multimaps \
                    --num-threads ${NUM_TOTAL_THREADS} \
                    ${SAMPLE_TSV} ${init_dir} ${SPECIES_PARA[@]}
    cd ${init_dir} && make sorted_reads | tee sargasso.log  2>&1
fi

####################################################################################
#### MAPPING
echo "Running Sargasso parameter tests...."
for mismatch in ${MISMATCH_SETTING}; do
    for minmatch in ${MINMATCH_SETTING}; do
        for multimap in ${MUTLIMAP_SETTING}; do
            echo "testing ${mismatch}_${minmatch}_${multimap}"
            out_dir=${OUTPUT_DIR}/${mismatch}_${minmatch}_${multimap}
            species_separator rnaseq --mapper-executable ${MAPPER_EXECUTABLE}  --sambamba-sort-tmp-dir=${TMP_DIR} \
                    --mismatch-threshold ${mismatch} --minmatch-threshold ${minmatch} --multimap-threshold ${multimap} --reject-multimaps \
                    --num-threads ${NUM_TOTAL_THREADS} \
                    ${SAMPLE_TSV} ${out_dir} ${SPECIES_PARA[@]}
            ## we link what we need
            ln -s -f ${init_dir}/mapper_indexes ${out_dir}
            ln -s -f ${init_dir}/raw_reads ${out_dir}
            ln -s -f ${init_dir}/mapped_reads ${out_dir}
            ln -s -f ${init_dir}/sorted_reads ${out_dir}
            cd ${out_dir} && make filtered_reads >sargasso.log 2>&1 &
        done
    done
    wait $(jobs -p)
done




echo "
library(dplyr)
library(ggplot2)
library(magrittr)
library(purrr)
library(readr)
library(reshape2)
library(stringr)

samples=c(\""`echo ${SAMPLES} | sed 's/ /","/g'`"\")
origin=c(\""`echo ${SAMPLES_ORIGIN} | sed 's/ /","/g'`"\") %>% set_names(samples)
result_dir=\""${OUTPUT_DIR}"\"
tb <- lapply(c("`echo ${MISMATCH_SETTING} | sed 's/ /,/g'`"),function(number_mismatch){
          lapply(c("`echo ${MINMATCH_SETTING} | sed 's/ /,/g'`"),function(min_match){
            lapply(c("`echo ${MUTLIMAP_SETTING} | sed 's/ /,/g'`"),function(number_multimap){
                file=file.path(result_dir,str_c(number_mismatch,'_',min_match,'_',number_multimap),'filtered_reads','overall_filtering_summary.txt')
                read_csv(file) %>% dplyr::select(Sample,contains('Reads')) %>%
                  mutate(Parameters=str_c(number_mismatch,min_match,number_multimap,sep = '_')) %>%
                  tidyr::pivot_longer(cols=contains('Reads'),names_to='type',values_to = 'count') %>%
                  mutate(origin=origin[Sample])
            }) %>% purrr::reduce(rbind)
          }) %>% purrr::reduce(rbind)
        }) %>% purrr::reduce(rbind)


count_table <- lapply(c(\""`echo ${SAMPLES} | sed 's/ /","/g'`"\")%>%set_names(.),function(sample){
  file<-list.files(path = \""${OUTPUT_DIR}"\", pattern = str_c(sample,'."${SAMPLES_ORIGIN%% *}".log.out',sep=''), all.files = FALSE,full.names = TRUE, recursive = TRUE,ignore.case = TRUE)[1]
  readLines(file,n= grep('Number of input reads',readLines(file))) %>% tail(1) %>% strsplit('\t') %>% extract2(1) %>% extract(2)
})

tb %<>% mutate(total_count=count_table[Sample] %>% unlist() %>% as.numeric()) %>%
  mutate(Percentage=count/total_count) %>%
  mutate_at(vars(one_of(c('Parameters','type'))),as.factor) %>%
  mutate(Parameters=factor(Parameters,level=unique(Parameters)))

.adjust_plot_size<-function(num_plots){
  num_row<-sqrt(num_plots) %>% ceiling()
  num_column<-num_plots/num_row %>% ceiling()
  c('width'= max(num_row/2,1),'height'=max(num_column/2,1))
}
sf <- .adjust_plot_size(length(origin)*3)

p_count <- tb %>%
  ggplot(aes(x=Parameters, y=count)) +
  geom_point(aes(color=Sample)) +
  facet_wrap(~origin + type,scales='free')+ coord_flip()
ggsave(file.path(result_dir,'counts.${PLOT_FORMAT}'),plot=p_count,width=12*sf,height=12*sf)

p_percentage <-tb %>%
  ggplot(aes(x=Parameters, y=Percentage)) +
  geom_point(aes(color=Sample)) +
  facet_wrap(~origin + type,scales='free')+ coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.01))
ggsave(file.path(result_dir,'percentage.${PLOT_FORMAT}'),plot=p_percentage,width=12*sf,height=12*sf)
" > ${OUTPUT_DIR}/plot.R

Rscript ${OUTPUT_DIR}/plot.R

echo "Test finish!"
echo "Please check result at: ${OUTPUT_DIR}/percentage.${PLOT_FORMAT} and ${OUTPUT_DIR}/counts.${PLOT_FORMAT}"



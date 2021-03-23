#!/usr/bin/env bash

set -o nounset
#set -o errexit
#set -o xtrace

function usage {
  cat <<EOT

sargasso_parameter_test.sh <data-type> <--samples-origin 'mouse mouse rat'>
[--mapper-executable STAR2.7.0f]
[--reads-base-dir /home/user/]]
[--mismatch-setting '0 2']
[--minmatch-setting '0 2']
[--multimap-setting '1']
[--num-threads 16]
[--plot-format png]
[--skip-init-run]
 <samples-file> <output-dir>
(<species> <species-info>)
(<species> <species-info>)
...

EOT
exit
}

subcommand=$1
case $subcommand in
    "" | "-h" | "--help") usage;;
    "rnaseq" |"dnaseq") DATA_TYPE=$1; shift;;
    * ) usage ;;
esac

OPTS="$(getopt -o h -l help,samples-origin:,skip-init-run,reads-base-dir:,mismatch-setting:,minmatch-setting:,multimap-setting:,mapper-executable:,num-threads:,plot-format: --name "$(basename "$0")" -- "$@")"
if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi
if [ "$#" -eq  "0"  ] ; then usage ; exit 1 ; fi

SKIP_INIT_RUN='no'
READS_BASE_DIR='/'
MISMATCH_SETTING='0 2 4 6 8 10'
MINMATCH_SETTING='0 2 5 10 '
MULTIMAP_SETTING='1'
MAPPER_EXECUTABLE=STAR
NUM_THREADS=1
PLOT_FORMAT=pdf
SPECIES_PARA=()

eval set -- "$OPTS"
while true; do
  case $1 in
    -h | --help ) usage;;
    --samples-origin ) SAMPLES_ORIGIN=$2; shift 2 ;;
    --skip-init-run ) SKIP_INIT_RUN='yes'; shift ;;
    --reads-base-dir ) READS_BASE_DIR=$2; shift 2 ;;
    --mismatch-setting ) MISMATCH_SETTING=$2; shift 2 ;;
    --minmatch-setting ) MINMATCH_SETTING=$2; shift 2 ;;
    --multimap-setting ) MULTIMAP_SETTING=$2; shift 2 ;;
    --mapper-executable ) MAPPER_EXECUTABLE=$2; shift 2 ;;
    --num-threads ) NUM_THREADS=$2; shift 2 ;;
    --plot-format ) PLOT_FORMAT=$2; shift 2 ;;
    -- ) shift; break  ;;
    * ) usage ;;
  esac
done
[[ -z ${SAMPLES_ORIGIN} ]] && echo "Error: --samples-origin is missing!" && exit 1

# We pass samples_file, output_dir and species_para
[[ "$#" -lt 6 ]] && usage && exit 1

SAMPLES_FILE=$1 && shift
[[ ! -s ${SAMPLES_FILE} ]] && echo "Error: sample_file(${SAMPLES_FILE}) is empty or cannot be found!" && exit 1
OUTPUT_DIR=$1 && shift
SPECIES_PARA+=($@)

## Extract samples from sample file
SAMPLES=`cut -d ' ' -f -1 ${SAMPLES_FILE} | paste -d " " -s`

# Check "samples" and "--samples-origin" have the same number of elements
number_sample=`echo "${SAMPLES}" | awk -F' ' '{print NF}'`
number_sample_origin=`echo "${SAMPLES_ORIGIN}" | awk -F' ' '{print NF}'`
[[ ${number_sample} -ne ${number_sample_origin} ]] && echo "Error: number of sample does not equal to number of sample origin." && exit 1

# Only support png and pdf for now
[[ ! "${PLOT_FORMAT}" =~ ^(pdf|png)$ ]] && echo "Error: --plot-format can only be one of pdf/png!" && exit 1


mkdir -p ${OUTPUT_DIR}

TMP_DIR=${OUTPUT_DIR}/tmp
init_dir=${OUTPUT_DIR}/init_run
mkdir -p ${TMP_DIR}


if [  "${SKIP_INIT_RUN}" == "no"  ]; then
    echo "Running Sargasso initial mapping...."
    # Execute the first sargasso run to get the mapped reads
    species_separator ${DATA_TYPE} --mapper-executable ${MAPPER_EXECUTABLE}  \
        --reads-base-dir=${READS_BASE_DIR} --sambamba-sort-tmp-dir=${TMP_DIR} \
        --mismatch-threshold 0 --minmatch-threshold 0 --multimap-threshold 1 --reject-multimaps \
        --num-threads ${NUM_THREADS} \
        ${SAMPLES_FILE} ${init_dir} ${SPECIES_PARA[@]}
    (cd ${init_dir} && make sorted_reads | tee sargasso.log  2>&1)
fi

echo "Running Sargasso parameter tests...."
for mismatch in ${MISMATCH_SETTING}; do
    for minmatch in ${MINMATCH_SETTING}; do
        for multimap in ${MULTIMAP_SETTING}; do
            echo "testing ${mismatch}_${minmatch}_${multimap}"
            out_dir=${OUTPUT_DIR}/${mismatch}_${minmatch}_${multimap}
            species_separator ${DATA_TYPE} --mapper-executable ${MAPPER_EXECUTABLE} \
               --reads-base-dir=${READS_BASE_DIR} --sambamba-sort-tmp-dir=${TMP_DIR} \
               --mismatch-threshold ${mismatch} --minmatch-threshold ${minmatch} --multimap-threshold ${multimap} --reject-multimaps \
               --num-threads ${NUM_THREADS} \
               ${SAMPLES_FILE} ${out_dir} ${SPECIES_PARA[@]}
            # Link what we need from the first Sargasso run
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
            lapply(c("`echo ${MULTIMAP_SETTING} | sed 's/ /,/g'`"),function(number_multimap){
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

echo "Making plots..."
Rscript ${OUTPUT_DIR}/plot.R  >${OUTPUT_DIR}/plot.log 2>&1

echo "Test finished!"
echo "Please check results: ${OUTPUT_DIR}/percentage.${PLOT_FORMAT} and ${OUTPUT_DIR}/counts.${PLOT_FORMAT}"


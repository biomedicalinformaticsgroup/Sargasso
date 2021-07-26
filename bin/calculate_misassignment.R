require(magrittr,quietly = T,warn.conflicts =F)
require(dplyr,quietly = T,warn.conflicts =F)
require(readr,quietly = T,warn.conflicts =F)
require(purrr,quietly = T,warn.conflicts =F)
require(stringr,quietly = T,warn.conflicts =F)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Rscript calculate_misassignment.R  'mouse' '01_HMR,02_HMR,03_HMR,04_HMR' 'human,mouse,rat' '09_Pc,10_Pc,11_Pc,12_Pc,05_EC,06_EC,07_EC,08_EC' 'rat;rat;rat;rat;human;human;human;human' 0 '/home/xinhe/Projects/bbb_pbaxter/misassignment/gene_length.csv' '/home/xinhe/Projects/bbb_pbaxter/results/sargasso/filtered_reads/overall_filtering_summary.txt' '/home/xinhe/Projects/bbb_pbaxter/results/read_counts/' '/home/xinhe/Projects/bbb_pbaxter/results/misassignment'")
}

species_of_interest=args[1]
target_samples=args[2]
target_species=args[3]
reference_samples=args[4]
reference_species=args[5]
paired=as.numeric(args[6])
gene_lengths_file<-file.path(args[7])
overall_filtering_summary_file=file.path(args[8])
read_count_dir=file.path(args[9])
result_dir=file.path(args[10])



# species_of_interest='mouse'  # the species we want to eastimate the misassignment
# target_samples="01_HMR,02_HMR,03_HMR,04_HMR"
# target_species="human,mouse,rat"
# reference_samples="09_Pc,10_Pc,11_Pc,12_Pc,05_EC,06_EC,07_EC,08_EC"
# reference_species="rat;rat;rat;rat;human;human;human;human"
# paired=FALSE
# gene_lengths_file=file.path('misassignment/gene_length.csv')
# overall_filtering_summary_file=file.path('results/sargasso/filtered_reads/overall_filtering_summary.txt')
# read_count_dir=file.path('results/read_counts/')
# result_dir = file.path('results/misassignment')

# we check parameter
if(!length(args)==10){
    stop('all parameters are required.')
}

n_ref_samples <- strsplit(reference_samples,',') %>% unlist() %>% length()
n_ref_species <- strsplit(reference_species,';') %>% unlist() %>% length()
if(! n_ref_samples == n_ref_species){
    stop('each reference samples should have a matching reference species.')
}

n_target_sample <- strsplit(target_samples,',') %>% unlist() %>% length()
if(paired && !n_target_sample == n_ref_samples){
    stop('each target sample should have a matching reference sample in the paired sample setup.')
}

if(!file.exists(gene_lengths_file)){
    stop(str_c('gene_length file does not exist.'))
}

if(!file.exists(overall_filtering_summary_file)){
    stop(str_c('sargasso overall_filtering_summary file does not exist.'))
}

if(!dir.exists(read_count_dir)){
    stop(str_c('read_count_dir does not exist.'))
}

if(!file.exists(result_dir)){
    dir.create(result_dir)
}

# read in gene lenth dataframe
gene_lengths<-read_csv(gene_lengths_file)






calculate_per_gene_misassignment_percentages <- function(species_of_interest,
target_samples,target_species,
reference_samples,reference_species,
paired=FALSE,
gene_lengths){

    target_samples %<>% strsplit(.,',') %>% unlist()
    target_species %<>% strsplit(.,',') %>% unlist()
    reference_samples  %<>% strsplit(.,',') %>% unlist()
    reference_species  %<>% strsplit(.,';') %>% unlist()

    # we workout all species involved
    all_species=c(target_species,
    species_of_interest,
    strsplit(reference_species,',') %>% unlist()
    ) %>% unique

    # For each 'species_of_interest' gene, calculate an approximate gene “expression” in 'reference_samples'(rat)
    # containing only mouse and human material, measured in fragments per kilobase
    # per million mapped reads: FPKM_mm_hs(i)= 10^9*c(i) / (l(i)*R)
    # where c(i) is the number of reads assigned to rat gene in the mouse plus human
    # samples, l(i) is the length of the longest transcript of the gene, and R is the
    # total number of reads assigned to all genes of all species in the mouse plus
    # human samples.

    ## step 30 in paper
    # we group the reference sample by reference species
    # $human
    # [1] "05_EC" "06_EC" "07_EC" "08_EC"
    #
    # $mouse
    # [1] "13_As"  "14_As" "15_As_replacement" "16_As"
    refs <- split(reference_samples,reference_species)

    ## If the reference samples are paired with the target samples,
    ## we check if they have the same number of sample.
    if(paired){
        sample_all_paired = lapply(refs,function(x){
            length(x)==length(target_samples)}
        ) %>%
            unlist() %>% all

        if(!sample_all_paired)
        stop('Difference number of reference samples / target samples, cannot match pair.')

        ## We print out the pairing
        lapply(1:length(target_samples),function(i){
            c(setNames(c(target_samples[i]),str_c(target_species,collapse = ',')),sapply(refs,extract2,i)) %>% t() %>% as.data.frame()
        })%>%reduce(rbind) %>% print()
    }

    # for each group of reference samples, we work out the species_of_interest fpkm (misassignment)
    species_of_interest_in_reference_sample_fpkm <- refs %>% lapply(function(rsa){ # rsa: reference_sample
        # gene               `05_EC` `06_EC` `07_EC` `08_EC` counts gene_length max_transcript_length
        # <chr>                <dbl>   <dbl>   <dbl>   <dbl>  <dbl>       <dbl>                 <dbl>
        # 1 ENSRNOG00000049070     315     522     630     663   2130        1025                  1025
        # 2 ENSRNOG00000054139     189     197     210     294    890        2325                  2325
        # 3 ENSRNOG00000054846      75      91     143     114    423        1527                  1527
        ref_gene_counts_and_lengths <- get_gene_counts_and_lengths(
        rsa, species_of_interest, gene_lengths, sum_counts = TRUE)

        # 05_EC     06_EC     07_EC     08_EC    counts
        # 36191427  31824279  41255399  46717964 155989069
        ref_total_reads <- calculate_total_reads_for_species(
        rsa, all_species, sum_counts = TRUE)

        # This is using the total number of reads(counts) to workout the fpkm.
        # The per-sample fpkm is just for debugging purpose
        # gene               `05_EC` `06_EC` `07_EC` `08_EC`  counts
        # <chr>                <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
        # 1 ENSRNOG00000049070  8.49   16.0    14.9    13.8    13.3
        # 2 ENSRNOG00000055837  3.55    3.36    3.71    3.01    3.39
        # 3 ENSRNOG00000054139  2.25    2.66    2.19    2.71    2.45
        # 4 ENSRNOG00000054846  1.36    1.87    2.27    1.60    1.78
        # 5 ENSRNOG00000046657  1.69    1.69    1.52    0.984   1.43
        # 6 ENSRNOG00000011175  0.130   0.148   0.0916  0.222   0.151
        # 7 ENSRNOG00000030807  0.0285  0.162   0.125   0.199   0.132
        c('counts',rsa) %>%
            calculate_fpkm( ref_gene_counts_and_lengths, ref_total_reads) %>%
            dplyr::arrange(gene)
    })

    ## step 31 on paper
    # work out the species_of_interest fpkm in the target samples
    # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR` counts gene_length max_transcript_length
    # <chr>                 <dbl>    <dbl>    <dbl>    <dbl>  <dbl>       <dbl>                 <dbl>
    # 1 ENSRNOG00000000001       52       41       44       38    175        1416                  1416
    # 2 ENSRNOG00000000007        0        0        1        0      1        3866                  3381
    # 3 ENSRNOG00000000008      593      726      808     1093   3220        1747                  1747
    # 4 ENSRNOG00000000009        1        4        2        0      7        1361                  1361
    # 5 ENSRNOG00000000010        4        0        0        0      4        2444                  2444
    target_gene_counts_and_lengths <- get_gene_counts_and_lengths(
    target_samples, species_of_interest, gene_lengths, sum_counts = TRUE) %>%
    dplyr::arrange(gene)

    # 01_HMR    02_HMR    03_HMR    04_HMR
    # 149055415 122723195  77948579 115554113
    target_total_reads <- calculate_total_reads_for_species(
    target_samples,species_of_interest, sum_counts = FALSE)

    # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`
    # <chr>                 <dbl>    <dbl>    <dbl>    <dbl>
    # 1 ENSRNOG00000000001   1.63      1.48   0.948      0.881
    # 2 ENSRNOG00000000007   0         0      0.00902    0
    # 3 ENSRNOG00000000008  15.0      21.2   14.1       20.5
    # 4 ENSRNOG00000000009   0.0325    0.150  0.0448     0
    # 5 ENSRNOG00000000010   0.0725    0      0          0
    # 6 ENSRNOG00000000012   0.178     0      0.215      0
    # 7 ENSRNOG00000000017  16.5      12.1   15.9       10.3
    species_of_interest_in_target_sample_fpkm <- target_samples %>%
        calculate_fpkm(target_gene_counts_and_lengths, target_total_reads) %>%
        arrange(gene)


    ## step 32 on paper
    ## For each target sample, calculate the approximate ratio of
    ## the reference_species ('mouse plus human’) to the species_of_interest (rat) RNA in the sample.
    sargasso_filtering_summary <-  overall_filtering_summary_file %>% read_csv


    ## for each reference species, we workout a P(% of misassignment)
    P<-names(species_of_interest_in_reference_sample_fpkm) %>% set_names(.) %>% lapply(function(s){
        # refrence species, this will hopefully handle the mono/co reference samples
        rsps <- strsplit(s,',') %>% unlist()
        rsa <- refs[[s]] # reference samlple

        #FPKM_mm_hsi, depend on 'paired' parameter, avg OR the per-sample fpkm will be used
        # gene                 avg `05_EC` `06_EC` `07_EC` `08_EC`
        # <chr>              <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
        # 1 ENSRNOG00000000001     0       0       0       0       0
        # 2 ENSRNOG00000000007     0       0       0       0       0
        # 3 ENSRNOG00000000008     0       0       0       0       0
        f1 <- species_of_interest_in_reference_sample_fpkm[[s]] %>%
        dplyr::rename(avg=counts)

        # d
        # 04_HMR    02_HMR    03_HMR    01_HMR
        # 0.4278663 0.2358971 0.8291291 0.2187826
        d <- sargasso_filtering_summary %>%
            dplyr::select(str_c('Assigned-Reads-',species_of_interest),str_c('Assigned-Reads-',rsps)) %>%
            mutate(toTargetRatio = (rowSums(.)-.[[1]])/.[[1]]) %>%
            dplyr::select(toTargetRatio) %>%
            mutate(Sample=sargasso_filtering_summary$Sample) %>%
            dplyr::select(Sample, toTargetRatio) %>%
            filter(Sample %in% target_samples) %>%
            tibble::deframe() %>% extract(target_samples)


        # FPKM_rni
        # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`
        # <chr>                 <dbl>    <dbl>    <dbl>    <dbl>
        # 1 ENSRNOG00000000001   1.63      1.48   0.948      0.881
        # 2 ENSRNOG00000000007   0         0      0.00902    0
        # 3 ENSRNOG00000000008  15.0      21.2   14.1       20.5
        # 4 ENSRNOG00000000009   0.0325    0.150  0.0448     0
        # 5 ENSRNOG00000000010   0.0725    0      0          0
        f2 <- species_of_interest_in_target_sample_fpkm

        if (!all(f1$gene == f2$gene)) {
            stop('gene in reference are not the same as gene in target.')
        }

        if(paired){
            ## for each refernce, we have 1 fpkm
            #                       05_EC 06_EC 07_EC 08_EC
            # ENSRNOG00000000001     0     0     0     0
            # ENSRNOG00000000007     0     0     0     0
            # ENSRNOG00000000008     0     0     0     0
            # ENSRNOG00000000009     0     0     0     0
            # ENSRNOG00000000010     0     0     0     0
            # ENSRNOG00000000012     0     0     0     0
            fpkm1 = f1 %>% dplyr::select(-avg) %>%
                tibble::column_to_rownames(var = 'gene') %>%
                as.matrix()
            numerator <- (( fpkm1 %>% as.matrix() %>% t()) * d )  %>% t()
        }else{
            ## we have one fpkm for all ref samples
            #                     avg
            # ENSRNOG00000000001   0
            # ENSRNOG00000000007   0
            # ENSRNOG00000000008   0
            # ENSRNOG00000000009   0
            # ENSRNOG00000000010   0
            # ENSRNOG00000000012   0
            fpkm1 = f1 %>% dplyr::select(gene,avg) %>%
                tibble::column_to_rownames(var = 'gene')
            numerator <- (fpkm1 %>% as.matrix()) %*% (d %>% as.matrix() %>% t())
        }


        #                         01_HMR     02_HMR       03_HMR     04_HMR
        # ENSRNOG00000000001  1.62657851  1.4761767  0.947681184  0.8813083
        # ENSRNOG00000000007  0.00000000  0.0000000  0.009020439  0.0000000
        # ENSRNOG00000000008 15.03476793 21.1866088 14.105591117 20.5463534
        # ENSRNOG00000000009  0.03254444  0.1498372  0.044817198  0.0000000
        # ENSRNOG00000000010  0.07249261  0.0000000  0.000000000  0.0000000
        # ENSRNOG00000000012  0.17842088  0.0000000  0.214991667  0.0000000
        denominator <- f2 %>%
            tibble::column_to_rownames(var = 'gene') %>%
            dplyr::select(names(d)) %>%
            as.matrix()

        # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`     p
        # <chr>                 <dbl>    <dbl>    <dbl>    <dbl> <dbl>
        # 1 ENSRNOG00000000001        0        0        0        0     0
        # 2 ENSRNOG00000000007      NaN      NaN        0      NaN     0
        # 3 ENSRNOG00000000008        0        0        0        0     0
        # 4 ENSRNOG00000000009        0        0        0      NaN     0
        # 5 ENSRNOG00000000010        0      NaN      NaN      NaN     0
        per_sample_p <- (100 * numerator /  denominator) %>%
            as.data.frame() %>% setNames(target_samples) %>%
            tibble::rownames_to_column('gene') %>%
            mutate_at(.vars = target_samples, list(~replace(., is.infinite(.), NaN))) %>%
            mutate(misprec = rowMeans(.[,target_samples], na.rm=TRUE)) %>%
            as_tibble()

        list(p=per_sample_p,d=d)
    })

    # $human
    # $human$p
    # # A tibble: 32,883 x 6
    # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`     p
    # <chr>                 <dbl>    <dbl>    <dbl>    <dbl> <dbl>
    #   1 ENSRNOG00000000001        0        0        0        0     0
    # 2 ENSRNOG00000000007      NaN      NaN        0      NaN     0
    # 3 ENSRNOG00000000008        0        0        0        0     0
    # 4 ENSRNOG00000000009        0        0        0      NaN     0
    # 5 ENSRNOG00000000010        0      NaN      NaN      NaN     0
    # # … with 32,873 more rows
    #
    # $human$d
    # 01_HMR     02_HMR     03_HMR     04_HMR
    # 0.45306509 0.45660661 0.06097529 0.13513955
    #
    # $mouse
    # $mouse$p
    # # A tibble: 32,883 x 6
    # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`     p
    # <chr>                 <dbl>    <dbl>    <dbl>    <dbl> <dbl>
    #   1 ENSRNOG00000000001        0        0        0        0     0
    # 2 ENSRNOG00000000007      NaN      NaN        0      NaN     0
    # 3 ENSRNOG00000000008        0        0        0        0     0
    # 4 ENSRNOG00000000009        0        0        0      NaN     0
    # 5 ENSRNOG00000000010        0      NaN      NaN      NaN     0
    # # … with 32,873 more rows
    #
    # $mouse$d
    # 01_HMR   02_HMR   03_HMR   04_HMR
    # 4.570748 4.239137 1.206085 2.337179

    # gene               p.human p.mouse     p
    # <chr>                <dbl>   <dbl> <dbl>
    # 1 ENSRNOG00000000001       0       0     0
    # 2 ENSRNOG00000000007       0       0     0
    # 3 ENSRNOG00000000008       0       0     0
    # 4 ENSRNOG00000000009       0       0     0
    # 5 ENSRNOG00000000010       0       0     0
    P_df<-P %>% lapply(extract2,'p') %>%
        purrr::reduce(left_join,by='gene',suffix=c(str_c('.',names(.)))) %>%
        dplyr::select(gene,starts_with('misprec')) %>%
        mutate(misprec=rowSums(.[2:ncol(.)]))

    # gene               human_counts `05_EC` `06_EC` `07_EC` `08_EC` mouse_counts `13_As` `14_As` `15_As_replacement` `16_As` `01_HMR` `02_HMR` `03_HMR` `04_HMR`
    # <chr>                     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>        <dbl>   <dbl>   <dbl>               <dbl>   <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
    # 1 ENSRNOG00000000001            0       0       0       0       0            0       0       0                   0       0   1.63      1.48   0.948      0.881
    # 2 ENSRNOG00000000007            0       0       0       0       0            0       0       0                   0       0   0         0      0.00902    0
    # 3 ENSRNOG00000000008            0       0       0       0       0            0       0       0                   0       0  15.0      21.2   14.1       20.5
    debug_df<-species_of_interest_in_reference_sample_fpkm %>% map2(names(.),.,function(s,df){
        df %>% dplyr::rename_at(vars(counts),funs(str_c(s,'_',.)))
    }) %>% purrr::reduce(left_join) %>% left_join(species_of_interest_in_target_sample_fpkm)

    # $human
    # 04_HMR     02_HMR     03_HMR     01_HMR
    # 0.13513955 0.45660661 0.06097529 0.45306509
    #
    # $mouse
    # 04_HMR   02_HMR   03_HMR   01_HMR
    # 2.337179 4.239137 1.206085 4.570748
    D_df <- P %>% lapply(extract2,'d')

    ## we work out the composition of the reference sample
    # refs_composition <- lapply(refs,function(ss){
    refs_composition <- sargasso_filtering_summary %>%
    # dplyr::filter(Sample %in% ss) %>%
        dplyr::select(Sample,starts_with('Assigned-Reads')) %>%
        tidyr::pivot_longer(-Sample,names_to='type') %>%
        group_by(Sample) %>%
        mutate(prec=round(100*value/sum(value),digits=1)) %>%
        tidyr::pivot_wider(Sample,names_from = type,values_from = prec ) %>%
        dplyr::rename_at(vars(-Sample),gsub,pattern='Assigned-Reads-',replacement='',fixed=TRUE) %>%
        ungroup()
    # })


    list(
    P=P_df %>% left_join(debug_df),
    D=D_df,
    refs_composition=refs_composition
    )
}


read_counts <- function(sample, species) {
    counts_file_name <- file.path(read_count_dir, str_c(sample, ".", species,".counts"))
    counts_file_name %>% read_tsv(col_names=c("gene", str_c(sample)))
}

get_gene_counts_and_lengths <- function(samples, species, gene_lengths, sum_counts = FALSE) {
    res <- samples %>%
        map(read_counts, species) %>%
        purrr::reduce(inner_join)

    if (sum_counts) {
        res %<>% mutate(counts = rowSums(.[,-1]))
    }

    res %>% inner_join(gene_lengths)
}

calculate_total_reads_for_species <- function(samples, species, sum_counts = FALSE) {
    res <- samples %>%
    map_dbl(function(sample) {
        species %>% map(~read_counts(sample, .) %>%
            dplyr::pull(2) %>%
            sum()) %>%
            purrr::reduce(sum)
    })

    if (sum_counts) {
        res %<>% set_names(samples)
        res['counts']=sum(res)
    } else {
        res %<>% set_names(samples)
    }
    res
}


calculate_fpkm <- function(samples, gene_counts_and_lengths, total_reads){
    samples %>%
        map(function(sa) {
            read_counts <- gene_counts_and_lengths %>% pull(sa)
            gene_lengths <- gene_counts_and_lengths %>% pull(max_transcript_length)
            gene_counts_and_lengths %>%
                dplyr::select(gene) %>%
                mutate(!!sa:=10^9 * as.double(read_counts) /
                    as.double(gene_lengths) /
                    as.double(total_reads[[sa]]))
        }) %>% reduce(inner_join)
}

calculate_species_ratios <- function(target_species){
    sargasso_filtering_summary <- read_overall_filtering_summary()

    tmp <- sargasso_filtering_summary %>% dplyr::select(starts_with('Assigned-Reads'))
    index <- match(str_c('Assigned-Reads-', target_species), names(tmp))

    tmp  %>%
        mutate(toTargetRatio = (rowSums(.)-.[[index]])/.[[index]]) %>%
        dplyr::select(toTargetRatio) %>%
        mutate(Sample=sargasso_filtering_summary$Sample) %>%
        dplyr::select(Sample, toTargetRatio)
}






res <- calculate_per_gene_misassignment_percentages(species_of_interest,
target_samples,target_species,
reference_samples,reference_species,
paired,
gene_lengths)

res$P %>% dplyr::select(gene,misprec,starts_with('misprec.')) %>% write.csv(file.path(result_dir,'gene_misassignment_percentage.csv'))
res$D %>% as.data.frame() %>% write.csv(file.path(result_dir,'sample_misassignment_percentage.csv'),)
res$refs_composition%>% write.csv(file.path(result_dir,'reference_sample_composition.csv'))
cat('Calculation finish. Results are stored in ')
cat(result_dir)

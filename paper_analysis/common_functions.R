library(DESeq2)
library(dplyr)
library(gplots)
library(lazyeval)
library(magrittr)
library(purrr)
library(RColorBrewer)
library(readr)
library(stringr)

rename_cols <- function(in_df, old_names, new_names) {
  out_df <- in_df %>% dplyr::rename_(.dots=setNames(old_names, new_names))
  return(out_df)
}

select_cols <- function(in_df, col_names) {
  out_df <- in_df %>% dplyr::select_(.dots=col_names)
  return(out_df)
}

remove_cols <- function(in_df, col_names) {
  out_df <- in_df %>% dplyr::select_(.dots=str_c("-", col_names))
  return(out_df)
}

quote_cols <- function(col_names) {
 return(str_c("`", col_names, "`")) 
}

get_gene_info <- function(species) {
    gene_info <- read_tsv(str_c("data/", species, "_ensembl-84/genes.tsv"), 
                          col_names = c("gene", "description", "chromosome", "gene_name"),
                          col_types = list(`chromosome` = col_character()))

    return(gene_info)
}

# Read in list of rat-mouse orthologous genes downloaded from Ensembl BioMart with the
# following query:
#
#  <Dataset name = "rnorvegicus_gene_ensembl" interface = "default" >
#  <Filter name = "with_homolog_mmus" excluded = "0"/>
#  <Attribute name = "ensembl_gene_id" />
#  <Attribute name = "mmusculus_homolog_ensembl_gene" />
#  <Attribute name = "mmusculus_homolog_orthology_type" />
#  </Dataset>
get_ortholog_info <- function() {
    orthologs <- read_tsv("data/rat_ensembl-84/mouse_orthologs.tsv",
                          col_names=c("rat_gene", "mouse_gene", "type"))

    remove_duplicate_genes <- function(gene_info, gene_column) {
      duplicates_removed <- gene_info %>%
        group_by_(gene_column) %>% 
        mutate(num_genes=n()) %>% 
        filter(num_genes == 1) %>%
        dplyr::select(-num_genes) %>%
        ungroup()
      
      return(duplicates_removed)
    }

    # Select only 1-to-1 orthologs
    o2o_orthologs <- orthologs %>%
      filter(type == "ortholog_one2one") %>%
      dplyr::select(-type) %>%
      mutate(o2o_ortholog="Y") %>%
      remove_duplicate_genes("mouse_gene") %>% 
      remove_duplicate_genes("rat_gene")

    mouse_genes <- get_gene_info("mouse")
    rat_genes <- get_gene_info("rat")

    same_name_genes <- (mouse_genes %>% 
                          dplyr::select(-chromosome, -description) %>% 
                          dplyr::rename(mouse_gene=gene)) %>% 
      inner_join(rat_genes %>% 
                   dplyr::select(-chromosome, -description) %>% 
                   dplyr::rename(rat_gene=gene)) %>% 
      dplyr::select(-gene_name) %>%
      mutate(same_name="Y") %>%
      remove_duplicate_genes("mouse_gene") %>% 
      remove_duplicate_genes("rat_gene")

    all_orthologs <- same_name_genes %>% merge(o2o_orthologs, all=T)
    all_orthologs[is.na(all_orthologs)] <- ""

    duplicated_mouse_genes <- (all_orthologs %>% 
                                 group_by(mouse_gene) %>% 
                                 mutate(num_genes=n()) %>% 
                                 filter(num_genes > 1) %>% 
                                 distinct())$mouse_gene

    all_orthologs %<>% filter(not(mouse_gene %in% duplicated_mouse_genes & 
                                  same_name == "Y"))

    duplicated_rat_genes <- (all_orthologs %>% 
                               group_by(rat_gene) %>% 
                               mutate(num_genes=n()) %>% 
                               filter(num_genes > 1) %>% 
                               distinct())$rat_gene

    all_orthologs %<>% filter(not(rat_gene %in% duplicated_rat_genes & 
                                  same_name == "Y"))

    return(all_orthologs)
}

read_counts <- function(counts_dir, sample, count_type, sample_suffix="") {
    counts_file_name <- str_c(counts_dir, "/", sample, ".", count_type, ".count")
    counts <- read_tsv(counts_file_name, col_names=c("gene", str_c(sample, sample_suffix)))
    return(counts)
}

get_samples <- function(conditions, replicates) {
  num_conditions <- length(conditions)
  num_replicates <- length(replicates)
  
  condition_vals <- c(sapply(conditions, rep, times=num_replicates))
  sample_vals <- rep(replicates, num_conditions)
  samples <- str_c(c(outer(replicates, conditions, str_c)))
  
  return(samples)
}

assemble_count_data <- function(counts_dir, conditions, replicates, count_type, 
                                sample_prefix="", sample_suffix="") {
  
  samples <- get_samples(conditions, replicates)
  count_data <- read_counts(counts_dir, samples[1], count_type, sample_suffix)
  samples <- samples[-1]
  
  for (sample in samples) {
    counts <- read_counts(counts_dir, sample, count_type, sample_suffix)
    count_data %<>% inner_join(counts)
  }
  
  row.names(count_data) <- count_data$gene
  count_data %<>% dplyr::select(-gene)

  return(count_data)
}

assemble_sample_data <- function(conditions, replicates, sample_prefix="") {
    num_conditions <- length(conditions)
    num_replicates <- length(replicates)

    condition_vals <- c(sapply(conditions, rep, times=num_replicates))
    sample_vals <- rep(replicates, num_conditions)
    row_names_vals <- str_c(sample_prefix, c(outer(replicates, conditions, str_c)))

    sample_data <- data.frame(condition=condition_vals,
                              sample=sample_vals,
                              row.names=row_names_vals)

    return(sample_data)
}

get_deseq2_dataset <- function(count_data, sample_data, filter_low_counts=TRUE, 
                               design_formula=~sample + condition) {
    
    dds <- DESeqDataSetFromMatrix(countData=count_data, colData=sample_data, design=design_formula)

    if (filter_low_counts) {
        dds <- dds[rowSums(counts(dds)) > 1, ]
    }

    return(DESeq(dds))
}

get_deseq2_results <- function(dds, condition, condition_base) {
    res <- results(dds, c("condition", condition, condition_base)) 
    
    print(summary(res))
    
    res %<>%
      as.data.frame() %>%
      tibble::rownames_to_column(var="gene") %>%
      dplyr::select(-baseMean, -lfcSE, -stat)

    return(res)
}

get_count_data <- function(dds, norm=T) {
  counts <- dds %>%
    counts(normalized=norm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="gene")

  return(counts)
}

plot_heat_map <- function(rld, sample_data) {
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    rownames(mat) <- colnames(mat) <- row.names(sample_data)
    hc <- hclust(distsRL)
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(mat, Rowv=as.dendrogram(hc), 
              symm=TRUE, trace="none",
              col = rev(hmcol), margin=c(5, 5))
}

get_fpkms <- function(all_counts, gene_lengths, samples, col_suffix) {
  all_counts %<>% inner_join(gene_lengths)
  
  for (sample in samples) {
    mmr <- sum(all_counts[[sample]]) / 1000000
    
    all_counts[,str_c(sample, col_suffix)] <- 
      map2_dbl(all_counts[[sample]], 
               all_counts[["max_transcript_length"]], 
               function(x, y) x / y / mmr * 1000)
  }
  
  all_counts %>% dplyr::select(gene, dplyr::contains(col_suffix))
}

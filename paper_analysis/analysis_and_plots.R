library(dplyr)
library(ggplot2)
library(magrittr)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(readr)
library(reshape2)
library(stringr)
library(topGO)

source("common_functions.R")

##### 

# FUNCTIONS

load_mapped_data <- function(csv_file, strategy) {
  data <- read_csv(csv_file) %>%
    rename_(.dots=setNames(
      list("correct_mapped_seq_count", "correct_mapped_seq_frac",
           "incorrect_mapped_seq_count", "incorrect_mapped_seq_frac"),
      c(str_c("correct_mapped_seq_count.", strategy),
        str_c("correct_mapped_seq_frac.", strategy),
        str_c("incorrect_mapped_seq_count.", strategy),
        str_c("incorrect_mapped_seq_frac.", strategy))))
  
  return(data)
}

load_assigned_data <- function(csv_file) {
  data <- read_delim(csv_file, " ",
                     col_names = c("gene", "conservative_assigned", 
                                   "best_assigned", "unfiltered_assigned"))
  
  data %<>% mutate(
    conservative_assigned_frac = conservative_assigned/unfiltered_assigned,
    best_assigned_frac = best_assigned/unfiltered_assigned)
  
  return(data)
}

collate_theoretical_results <- function(species) {
  base_dir <- "results/theoretical/"
  
  cons_results <- load_mapped_data(
    str_c(base_dir, "/", species, "_reads_per_gene_50.conservative.csv"), "conservative")
  best_results <- load_mapped_data(
    str_c(base_dir, "/",species, "_reads_per_gene_50.best.csv"), "best")

  assigned_results <- load_assigned_data(
    str_c(base_dir, "/assigned_", species, "_reads_per_gene_50.csv"))

  results <- cons_results %>% 
    inner_join(best_results) %>%
    inner_join(assigned_results) %>%
    left_join(get_gene_info(species)) %>%
    left_join(get_ortholog_info(), by=c("gene"=str_c(species, "_gene"))) %>% 
    dplyr::select(gene, gene_name, description, chromosome, 
           same_name, o2o_ortholog, seq_count,
           everything())
  
  return(results)
}

plot_mapped_reads_lost_histogram <- function(data, species) {
  hist_data <- data %>%
    dplyr::select(correct_mapped_seq_frac.conservative, correct_mapped_seq_frac.best) %>%
    dplyr::rename(Conservative=correct_mapped_seq_frac.conservative, Best=correct_mapped_seq_frac.best) %>%
    melt(id.vars=c(), variable.name="Strategy")
  
  ggplot(hist_data, aes(x=1-value, fill=Strategy)) + 
    geom_histogram(bins=21, boundary=0, position="dodge") +
    xlab("Fraction of mapped paired-end reads lost due to species separation") + 
    ylab("Number of genes") + 
    ggtitle(str_c("Theoretical reads covering", species, "transcriptome", sep=" "))
}

plot_mapped_reads_incorrect_histogram <- function(data, species) {
  hist_data <- data %>%
    dplyr::select(incorrect_mapped_seq_frac.conservative, incorrect_mapped_seq_frac.best) %>%
    dplyr::rename(Conservative=incorrect_mapped_seq_frac.conservative, Best=incorrect_mapped_seq_frac.best) %>%
    melt(id.vars=c(), variable.name="Strategy")
  
  ggplot(hist_data, aes(x=value, fill=Strategy)) + 
    geom_histogram(bins=21, boundary=0, position="dodge") +
    scale_x_log10() + 
    xlab("Fraction of mapped paired-end reads incorrectly assigned due to species separation") + 
    ylab("Number of genes") + 
    ggtitle(str_c("Theoretical reads covering", species, "transcriptome", sep=" "))
}

plot_assigned_reads_lost_histogram <- function(data, species) {
  hist_data <- data %>%
    dplyr::select(conservative_assigned_frac, best_assigned_frac) %>%
    dplyr::rename(Conservative=conservative_assigned_frac, Best=best_assigned_frac) %>%
    melt(id.vars=c(), variable.name="Strategy")
  
  ggplot(hist_data, aes(x=1-value, fill=Strategy)) + 
    geom_histogram(bins=21, boundary=0, position="dodge") +
    xlab("Fraction of feature assigned paired-end reads lost due to species separation") + 
    ylab("Number of genes") + 
    ggtitle(str_c("Theoretical reads covering", species, "transcriptome", sep=" "))
}

get_significant_genes <- function(term, GOdata) {
  genes_for_term <- genesInTerm(GOdata, term)[[1]]
  significant_genes <- sigGenes(GOdata)
  significant_genes_for_term <- intersect(genes_for_term, significant_genes)
  return(paste(sort(significant_genes_for_term), collapse=', '))
}

perform_go_analysis <- function(results, significant_genes, id_mapping) {
  gene_list <- factor(as.integer(results$gene %in% significant_genes$gene))
  names(gene_list) <- results$gene
  
  go_data <- new("topGOdata", ontology="BP", allGenes=gene_list,
                 annot=annFUN.org, mapping=id_mapping, ID="Ensembl")
  
  result_fisher <- runTest(go_data, algorithm="weight01", statistic="fisher")
  print(result_fisher)
  
  go_results <- GenTable(go_data, weight_fisher=result_fisher, orderBy="weight_fisher", topNodes=20)
  
  go_results$Genes <- sapply(go_results[,c('GO.ID')], function(x) get_significant_genes(x, go_data))
  
  return(go_results)
}

#####

mouse_results <- collate_theoretical_results("mouse")
rat_results <- collate_theoretical_results("rat")

plot_mapped_reads_lost_histogram(mouse_results, "mouse")
plot_mapped_reads_lost_histogram(rat_results, "rat")

mouse_results %>% nrow
mouse_results %>% filter(correct_mapped_seq_frac.conservative > 0.9) %>% nrow
mouse_results %>% filter(correct_mapped_seq_frac.best > 0.9) %>% nrow

rat_results %>% nrow
rat_results %>% filter(correct_mapped_seq_frac.conservative > 0.9) %>% nrow
rat_results %>% filter(correct_mapped_seq_frac.best > 0.9) %>% nrow

plot_mapped_reads_incorrect_histogram(mouse_results, "mouse")
plot_mapped_reads_incorrect_histogram(rat_results, "rat")

plot_assigned_reads_lost_histogram(mouse_results, "mouse")
plot_assigned_reads_lost_histogram(rat_results, "rat")

mouse_results %>% nrow
mouse_results %>% filter(conservative_assigned_frac < 0.9) %>% nrow
mouse_results %>% filter(conservative_assigned_frac < 0.5) %>% nrow

mouse_results %>% filter(best_assigned_frac < 0.9) %>% nrow
mouse_results %>% filter(best_assigned_frac < 0.5) %>% nrow

mouse_go_results.conservative <- perform_go_analysis(
  mouse_results, mouse_results %>% filter(conservative_assigned_frac < 0.8), "org.Mm.eg.db") 
mouse_go_results.best <- perform_go_analysis(
  mouse_results, mouse_results %>% filter(best_assigned_frac < 0.8), "org.Mm.eg.db")

rat_results %>% nrow
rat_results %>% filter(conservative_assigned_frac < 0.9) %>% nrow
rat_results %>% filter(conservative_assigned_frac < 0.5) %>% nrow

rat_results %>% filter(best_assigned_frac < 0.9) %>% nrow
rat_results %>% filter(best_assigned_frac < 0.5) %>% nrow

rat_go_results.conservative <- perform_go_analysis(
  rat_results, rat_results %>% filter(conservative_assigned_frac < 0.8), "org.Rn.eg.db") 
rat_go_results.best <- perform_go_analysis(
  rat_results, rat_results %>% filter(best_assigned_frac < 0.8), "org.Rn.eg.db")

mouse_results %>% filter(incorrect_mapped_seq_frac.conservative > 0.01) %>% nrow
mouse_results %>% filter(incorrect_mapped_seq_frac.best > 0.01) %>% nrow

mouse_go_results.incorrect <- perform_go_analysis(
  mouse_results, mouse_results %>% filter(incorrect_mapped_seq_frac.best > 0.01), "org.Mm.eg.db") 

rat_results %>% filter(incorrect_mapped_seq_frac.conservative > 0.01) %>% nrow
rat_results %>% filter(incorrect_mapped_seq_frac.best > 0.01) %>% nrow

rat_go_results.incorrect <- perform_go_analysis(
  rat_results, rat_results %>% filter(incorrect_mapped_seq_frac.best > 0.01), "org.Rn.eg.db") 

library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(stringr)

source("common_functions.R")

##### 

# FUNCTIONS

plot_reads_lost_histogram <- function(data, species) {
  hist_data <- data %>%
    select(correct_mapped_seq_frac.conservative, correct_mapped_seq_frac.best) %>%
    dplyr::rename(Conservative=correct_mapped_seq_frac.conservative, Best=correct_mapped_seq_frac.best) %>%
    melt(id.vars=c(), variable.name="Strategy")
  
  ggplot(hist_data, aes(x=1-value, fill=Strategy)) + 
    geom_histogram(bins=21, boundary=0, position="dodge") +
    xlab("% of paired-end reads lost due to species separation") + 
    ylab("Number of genes") + 
    ggtitle(str_c("Theoretical reads covering", species, "transcriptome", sep=" "))
}

load_theoretical_data <- function(csv_file, strategy) {
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

collate_theoretical_results <- function(cons_results, best_results, species) {
  results <- cons_results %>% 
    inner_join(best_results) %>%
    left_join(get_gene_info(species)) %>%
    left_join(get_ortholog_info(), by=c("gene"=str_c(species, "_gene"))) %>% 
    select(gene, gene_name, description, chromosome, 
           same_name, o2o_ortholog, seq_count,
           everything())
  
  return(results)
}

#####

mouse.conservative <- load_theoretical_data(
  "results/theoretical/mouse_reads_per_gene_50.conservative.csv", "conservative")
mouse.best <- load_theoretical_data(
  "results/theoretical/mouse_reads_per_gene_50.best.csv", "best")

rat.conservative <- load_theoretical_data(
  "results/theoretical/rat_reads_per_gene_50.conservative.csv", "conservative")
rat.best <- load_theoretical_data(
  "results/theoretical/rat_reads_per_gene_50.best.csv", "best")

mouse_results <- collate_theoretical_results(mouse.conservative, mouse.best, "mouse")
rat_results <- collate_theoretical_results(rat.conservative, rat.best, "rat")

plot_reads_lost_histogram(mouse_results, "mouse")
plot_reads_lost_histogram(rat_results, "rat")

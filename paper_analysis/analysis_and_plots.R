library(ggplot2)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(reshape2)
library(Rgraphviz)
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
    ggtitle(str_c("Theoretical reads covering", species, "transcriptome", sep=" ")) +
    theme_grey(base_size = 24) 
}

plot_mapped_reads_incorrect_histogram <- function(data, species, breaks) {
  hist_data <- data %>%
    dplyr::select(incorrect_mapped_seq_frac.conservative, incorrect_mapped_seq_frac.best) %>%
    dplyr::rename(Conservative=incorrect_mapped_seq_frac.conservative, Best=incorrect_mapped_seq_frac.best) %>%
    melt(id.vars=c(), variable.name="Strategy")
  
  ggplot(hist_data, aes(x=value, fill=Strategy)) + 
    geom_histogram(bins=21, boundary=0, position="dodge") +
    scale_x_log10(breaks=breaks, labels=as.character(breaks)) + 
    xlab("Fraction of mapped paired-end reads incorrectly assigned due to species separation") + 
    ylab("Number of genes") + 
    ggtitle(str_c("Theoretical reads covering", species, "transcriptome", sep=" ")) +
    theme_grey(base_size = 24) 
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
    ggtitle(str_c("Theoretical reads covering", species, "transcriptome", sep=" ")) +
    theme_grey(base_size = 24) 
}

get_significant_genes <- function(term, GOdata) {
  genes_for_term <- genesInTerm(GOdata, term)[[1]]
  significant_genes <- sigGenes(GOdata)
  significant_genes_for_term <- intersect(genes_for_term, significant_genes)
  return(significant_genes_for_term)
}

perform_go_analysis <- function(results, significant_genes, id_mapping) {
  gene_list <- factor(as.integer(results$gene %in% significant_genes$gene))
  names(gene_list) <- results$gene
  
  go_data <- new("topGOdata", ontology="BP", allGenes=gene_list,
                 annot=annFUN.org, mapping=id_mapping, ID="Ensembl")
  
  result_fisher <- runTest(go_data, algorithm="weight01", statistic="fisher")
  print(result_fisher)
  
  go_results <- GenTable(go_data, weight_fisher=result_fisher, orderBy="weight_fisher", topNodes=150)
  
  go_results$Genes <- sapply(go_results[,c('GO.ID')], function(x) get_significant_genes(x, go_data))
  
  return(list(go_data, result_fisher, go_results))
}

#####

mouse_results <- collate_theoretical_results("mouse")

100 * sum(mouse_results$correct_mapped_seq_count.conservative) / sum(mouse_results$seq_count)
100 * sum(mouse_results$correct_mapped_seq_count.best) / sum(mouse_results$seq_count)

rat_results <- collate_theoretical_results("rat")

100 * sum(rat_results$correct_mapped_seq_count.conservative) / sum(rat_results$seq_count)
100 * sum(rat_results$correct_mapped_seq_count.best) / sum(rat_results$seq_count)

plot_mapped_reads_lost_histogram(mouse_results, "mouse")
plot_mapped_reads_lost_histogram(rat_results, "rat")

100 * (mouse_results %>% filter(correct_mapped_seq_frac.conservative > 0.8) %>% nrow) / (mouse_results %>% nrow)
100 * (mouse_results %>% filter(correct_mapped_seq_frac.best > 0.8) %>% nrow) / (mouse_results %>% nrow)

100 * (rat_results %>% filter(correct_mapped_seq_frac.conservative > 0.8) %>% nrow) / (rat_results %>% nrow)
100 * (rat_results %>% filter(correct_mapped_seq_frac.best > 0.8) %>% nrow) / (rat_results %>% nrow)

plot_assigned_reads_lost_histogram(mouse_results, "mouse")
plot_assigned_reads_lost_histogram(rat_results, "rat")

100 * (mouse_results %>% filter(conservative_assigned_frac < 0.8) %>% nrow) / (mouse_results %>% nrow)
100 * (mouse_results %>% filter(best_assigned_frac < 0.8) %>% nrow) / (mouse_results %>% nrow)

100 * (rat_results %>% filter(conservative_assigned_frac < 0.8) %>% nrow) / (rat_results %>% nrow)
100 * (rat_results %>% filter(best_assigned_frac < 0.8) %>% nrow) / (rat_results %>% nrow)

mouse_go_results.conservative <- perform_go_analysis(
  mouse_results, mouse_results %>% filter(conservative_assigned_frac < 0.8), "org.Mm.eg.db") 
mouse_go_results.best <- perform_go_analysis(
  mouse_results, mouse_results %>% filter(best_assigned_frac < 0.8), "org.Mm.eg.db")

# GO term with most significant genes is "G-protein coupled receptor signaling pathway"
gpcr_genes <- mouse_go_results.conservative[[3]][2, "Genes"][[1]]

orthologs <- read_tsv("data/rat_ensembl-84/mouse_orthologs.tsv",
                      col_names=c("rat_gene", "mouse_gene", "type"))

# mouse genes with one-to-one orthologs
one_to_one <- orthologs %>% filter(type == "ortholog_one2one") %>% dplyr::select(mouse_gene) %>% distinct

# mouse genes with many-to-many orthologs
many_to_many <- orthologs %>% filter(type == "ortholog_many2many") %>% dplyr::select(mouse_gene) %>% distinct

# mouse genes with one mouse-many rat orthologs
one_mouse_to_many_rat <- orthologs %>% filter(type == "ortholog_one2many") %>% 
  group_by(mouse_gene) %>% mutate(num_genes=n()) %>% filter(num_genes > 1) %>% 
  dplyr::select(mouse_gene) %>% distinct

# number of GPCR genes with more than one rat ortholog
sum(gpcr_genes %in% many_to_many[[1]]) + sum(gpcr_genes %in% one_mouse_to_many_rat[[1]]) # 201
100 * (sum(gpcr_genes %in% many_to_many[[1]]) + sum(gpcr_genes %in% one_mouse_to_many_rat[[1]])) / length(gpcr_genes)
  
# number of all "conservative" genes (>20% loss) with more than one rat ortholog
all_cons_mouse_genes <- mouse_results %>% filter(conservative_assigned_frac < 0.8) %>% dplyr::select(gene) %>% extract2(1)
sum(all_cons_mouse_genes %in% many_to_many[[1]]) + sum(all_cons_mouse_genes %in% one_mouse_to_many_rat[[1]]) # 622 / 1153
100 * (sum(all_cons_mouse_genes %in% many_to_many[[1]]) + sum(all_cons_mouse_genes %in% one_mouse_to_many_rat[[1]])) / length(all_cons_mouse_genes)

# number of all "best" genes (>20% loss) with more than one mouse ortholog
all_best_mouse_genes <- mouse_results %>% filter(best_assigned_frac < 0.8) %>% dplyr::select(gene) %>% extract2(1)
sum(all_best_mouse_genes %in% many_to_many[[1]]) + sum(all_best_mouse_genes %in% one_mouse_to_many_rat[[1]]) # 4 / 88

# number of all "best" genes (>20% loss) with a one-to-one ortholog
sum(all_best_mouse_genes %in% one_to_one[[1]]) # 65 / 88
100 * sum(all_best_mouse_genes %in% one_to_one[[1]]) / length(all_best_mouse_genes)
  
#rat_go_results.conservative <- perform_go_analysis(
#  rat_results, rat_results %>% filter(conservative_assigned_frac < 0.8), "org.Rn.eg.db") 
#rat_go_results.best <- perform_go_analysis(
#  rat_results, rat_results %>% filter(best_assigned_frac < 0.8), "org.Rn.eg.db")

# GO term with most significant genes is "G-protein coupled receptor signaling pathway"
#rat_gpcr_genes <- rat_go_results.conservative[[3]][2, "Genes"][[1]]

# rat genes with one-to-one orthologs
#one_to_one_rat <- orthologs %>% filter(type == "ortholog_one2one") %>% dplyr::select(rat_gene) %>% distinct

# rat genes with many-to-many orthologs
many_to_many_rat <- orthologs %>% filter(type == "ortholog_many2many") %>% dplyr::select(rat_gene) %>% distinct

# rat genes with one rat-many mouse orthologs
#one_rat_to_many_mouse <- orthologs %>% filter(type == "ortholog_one2many") %>% 
#  group_by(rat_gene) %>% mutate(num_genes=n()) %>% filter(num_genes > 1) %>% 
#  dplyr::select(rat_gene) %>% distinct

# rat genes with one mouse-many rat orthologs
one_mouse_to_many_rat_RAT <- orthologs %>% filter(type == "ortholog_one2many") %>% 
  group_by(mouse_gene) %>% mutate(num_genes=n()) %>% filter(num_genes > 1) %>% ungroup %>%
  dplyr::select(rat_gene) %>% distinct

# number of all "conservative" genes (>20% loss) with more than one mouse ortholog
#all_cons_rat_genes <- rat_results %>% filter(conservative_assigned_frac < 0.8) %>% dplyr::select(gene) %>% extract2(1)
#sum(all_cons_rat_genes %in% many_to_many_rat[[1]]) + sum(all_cons_rat_genes %in% one_rat_to_many_mouse[[1]]) # 321 / 1212
#100 * (sum(all_cons_rat_genes %in% many_to_many_rat[[1]]) + sum(all_cons_rat_genes %in% one_rat_to_many_mouse[[1]])) / length(all_cons_rat_genes)

100 * sum(mouse_results$incorrect_mapped_seq_count.conservative) / sum(mouse_results$seq_count)
100 * sum(mouse_results$incorrect_mapped_seq_count.best) / sum(mouse_results$seq_count)

mouse_results %>% filter(incorrect_mapped_seq_frac.conservative > 0.01) %>% nrow
100 * (mouse_results %>% filter(incorrect_mapped_seq_frac.conservative > 0.01) %>% nrow) / (mouse_results %>% nrow)
mouse_results %>% filter(incorrect_mapped_seq_frac.best > 0.01) %>% nrow
100 * (mouse_results %>% filter(incorrect_mapped_seq_frac.best > 0.01) %>% nrow) / (mouse_results %>% nrow)

mouse_go_results.incorrect <- perform_go_analysis(
  mouse_results, mouse_results %>% filter(incorrect_mapped_seq_frac.best > 0.01), "org.Mm.eg.db") 

100 * sum(rat_results$incorrect_mapped_seq_count.conservative) / sum(rat_results$seq_count)
100 * sum(rat_results$incorrect_mapped_seq_count.best) / sum(rat_results$seq_count)

rat_results %>% filter(incorrect_mapped_seq_frac.conservative > 0.01) %>% nrow
100 * (rat_results %>% filter(incorrect_mapped_seq_frac.conservative > 0.01) %>% nrow) / (rat_results %>% nrow)
rat_results %>% filter(incorrect_mapped_seq_frac.best > 0.01) %>% nrow
100 * (rat_results %>% filter(incorrect_mapped_seq_frac.best > 0.01) %>% nrow) / (rat_results %>% nrow)
rat_results %>% filter(incorrect_mapped_seq_frac.best > 0.2) %>% nrow
100 * (rat_results %>% filter(incorrect_mapped_seq_frac.best > 0.2) %>% nrow) / (rat_results %>% nrow)

rat_go_results.incorrect <- perform_go_analysis(
  rat_results, rat_results %>% filter(incorrect_mapped_seq_frac.best > 0.01), "org.Rn.eg.db") 

all_best_incorrect_rat_genes <- rat_results %>% filter(incorrect_mapped_seq_frac.best > 0.01) %>% dplyr::select(gene) %>% extract2(1)
sum(all_best_incorrect_rat_genes %in% many_to_many_rat[[1]]) + sum(all_best_incorrect_rat_genes %in% one_mouse_to_many_rat_RAT[[1]]) # 321 / 1212
100 * (sum(all_best_incorrect_rat_genes %in% one_mouse_to_many_rat_RAT[[1]])) / length(all_best_incorrect_rat_genes)

plot_mapped_reads_incorrect_histogram(mouse_results, "mouse", c(0.0001, 0.001, 0.01, 0.1, 0.5))
plot_mapped_reads_incorrect_histogram(rat_results, "rat", c(0.0001, 0.001, 0.01, 0.1, 0.5, 1))


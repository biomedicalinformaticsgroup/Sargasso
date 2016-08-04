library(gridExtra)
library(VennDiagram)

source("common_functions.R")

##### 

# FUNCTIONS

get_up_reg_genes <- function(results) {
  genes <- results %>% 
    filter(padj < 0.1 & log2FoldChange > 0) %>% 
    extract2(1)
  
  return(genes)
}

get_down_reg_genes <- function(results) {
  genes <- results %>% 
    filter(padj < 0.1 & log2FoldChange < 0) %>% 
    extract2(1)
  
  return(genes)
}

draw_venn_diagram <- function(normal_genes, best_genes, conservative_genes) {
  normal_and_best <- intersect(normal_genes, best_genes)
  normal_and_conservative <- intersect(normal_genes, conservative_genes)
  best_and_conservative <- intersect(best_genes, conservative_genes)
  all <- intersect(normal_and_best, conservative_genes)
  
  venn_diag <- draw.triple.venn(
    length(normal_genes), length(best_genes), length(conservative_genes), 
    length(normal_and_best), length(best_and_conservative), 
    length(normal_and_conservative), length(all), 
    category=c("Normal", "Best", "Conservative"), 
    fill=c("red", "blue", "green"), alpha=0.5, ind=F,
    fontfamily = rep("helvetica", 7),
    cat.fontfamily=rep("helvetica", 3),
    cex=rep(3, 7), cat.cex=rep(3, 3))
  
  frame()
  grid.draw(venn_diag)
}

#####

conditions <- c("1", "2")
replicates <- c("AN", "BN", "CN")
counts_dir <- "results/mouse_activity/read_counts/"

# Do differential expression for the normal STAR mapping data
sample_data <- assemble_sample_data(conditions, replicates)
normal_count_data <- assemble_count_data(counts_dir, conditions, replicates, "normal")

normal_dds <- get_deseq2_dataset(normal_count_data, sample_data)
normal_res_2v1 <- get_deseq2_results(normal_dds, "2", "1")

normal_dds_all <- get_deseq2_dataset(
  normal_count_data, sample_data, filter_low_counts = FALSE)

# Do differential expression for the "conservative" SSS data
conservative_count_data <- assemble_count_data(counts_dir, conditions, replicates, "conservative")

conservative_dds <- get_deseq2_dataset(conservative_count_data, sample_data)
conservative_res_2v1 <- get_deseq2_results(conservative_dds, "2", "1")

conservative_dds_all <- get_deseq2_dataset(
  conservative_count_data, sample_data, filter_low_counts = FALSE)

# Do differential expression for the "best" SSS data
best_count_data <- assemble_count_data(counts_dir, conditions, replicates, "best")

best_dds <- get_deseq2_dataset(best_count_data, sample_data)
best_res_2v1 <- get_deseq2_results(best_dds, "2", "1")

best_dds_all <- get_deseq2_dataset(
  best_count_data, sample_data, filter_low_counts = FALSE)

#####

normal_up <- get_up_reg_genes(normal_res_2v1)
best_up <- get_up_reg_genes(best_res_2v1)
conservative_up <- get_up_reg_genes(conservative_res_2v1)

draw_venn_diagram(normal_up,best_up, conservative_up)
draw_venn_diagram(get_down_reg_genes(normal_res_2v1), 
                  get_down_reg_genes(best_res_2v1), 
                  get_down_reg_genes(conservative_res_2v1))

length(setdiff(normal_up, best_up))
length(setdiff(best_up, normal_up))

length(setdiff(normal_up, conservative_up))
length(setdiff(conservative_up, normal_up))

#####

counts_all <- get_count_data(normal_dds_all)

counts_best_all <- get_count_data(best_dds_all) %>%
  dplyr::rename(AN1_best=AN1, BN1_best=BN1, CN1_best=CN1,
                AN2_best=AN2, BN2_best=BN2, CN2_best=CN2)

counts_conservative_all <- get_count_data(conservative_dds_all) %>%
  dplyr::rename(AN1_conservative=AN1, BN1_conservative=BN1, CN1_conservative=CN1,
                AN2_conservative=AN2, BN2_conservative=BN2, CN2_conservative=CN2)

counts_all %<>% 
  inner_join(counts_best_all) %>%
  inner_join(counts_conservative_all)

#####

counts <- get_count_data(normal_dds)

counts_best <- get_count_data(best_dds) %>%
  dplyr::rename(AN1_best=AN1, BN1_best=BN1, CN1_best=CN1,
                AN2_best=AN2, BN2_best=BN2, CN2_best=CN2)

counts_conservative <- get_count_data(conservative_dds) %>%
  dplyr::rename(AN1_conservative=AN1, BN1_conservative=BN1, CN1_conservative=CN1,
                AN2_conservative=AN2, BN2_conservative=BN2, CN2_conservative=CN2)

results <- counts %>% 
    inner_join(counts_best) %>%
    inner_join(counts_conservative)

gene_lengths <- read_csv("results/gene_lengths.csv")

samples <- row.names(sample_data)

results %<>% left_join(
  counts_all %>% get_fpkms(gene_lengths, samples, "_fpkm"))

# Add differential expression results columns for normal STAR mapping
results %<>% mutate(l2fc_raw_normal = log2((AN2 + BN2 + CN2)/(AN1 + BN1 + CN1)))
results %<>% inner_join(normal_res_2v1, by="gene") %>%
    dplyr::rename(l2fc_normal = log2FoldChange,
                  pval_normal = pvalue,
                  padj_normal = padj)

# Add the same differential expression results columns for "best" strategy
results %<>% mutate(l2fc_raw_best = log2((AN2_best + BN2_best + CN2_best)/
                                             (AN1_best + BN1_best + CN1_best)))
results %<>% inner_join(best_res_2v1, by="gene") %>%
  dplyr::rename(l2fc_best = log2FoldChange,
                pval_best = pvalue,
                padj_best = padj)

# Add the same differential expression results columns for "conservative" strategy
results %<>% mutate(l2fc_raw_conservative = log2((AN2_conservative + BN2_conservative + CN2_conservative)/
                                           (AN1_conservative + BN1_conservative + CN1_conservative)))
results %<>% inner_join(conservative_res_2v1, by="gene") %>%
  dplyr::rename(l2fc_conservative = log2FoldChange,
                pval_conservative = pvalue,
                padj_conservative = padj)

# Add gene and ortholog info to results object
results %<>% left_join(get_gene_info("mouse")) %>%
  left_join(get_ortholog_info(), by=c("gene" = "mouse_gene"))

# Add gene and maximum transcript length information
results %<>% inner_join(gene_lengths)

#####

total_count_data <-
  read_counts(counts_dir, "AN1", "normal") %>%
  inner_join(read_counts(counts_dir, "BN1", "normal")) %>% 
  inner_join(read_counts(counts_dir, "CN1", "normal")) %>%
  inner_join(read_counts(counts_dir, "AN2", "normal")) %>%
  inner_join(read_counts(counts_dir, "BN2", "normal")) %>%
  inner_join(read_counts(counts_dir, "CN2", "normal")) %>%
  inner_join(read_counts(counts_dir, "AN1", "best", "_best")) %>%
  inner_join(read_counts(counts_dir, "BN1", "best", "_best")) %>%
  inner_join(read_counts(counts_dir, "CN1", "best", "_best")) %>%
  inner_join(read_counts(counts_dir, "AN2", "best", "_best")) %>%
  inner_join(read_counts(counts_dir, "BN2", "best", "_best")) %>%
  inner_join(read_counts(counts_dir, "CN2", "best", "_best")) %>%
  inner_join(read_counts(counts_dir, "AN1", "conservative", "_conservative")) %>%
  inner_join(read_counts(counts_dir, "BN1", "conservative", "_conservative")) %>%
  inner_join(read_counts(counts_dir, "CN1", "conservative", "_conservative")) %>%
  inner_join(read_counts(counts_dir, "AN2", "conservative", "_conservative")) %>%
  inner_join(read_counts(counts_dir, "BN2", "conservative", "_conservative")) %>%
  inner_join(read_counts(counts_dir, "CN2", "conservative", "_conservative"))

row.names(total_count_data) <- total_count_data$gene
total_count_data %<>% dplyr::select(-gene)

total_sample_data <- assemble_sample_data(
    c("1", "2", "1_best", "2_best", "1_conservative", "2_conservative"), c("AN", "BN", "CN"))

total_dds <- get_deseq2_dataset(total_count_data, total_sample_data)
total_rld <- rlog(total_dds)

plot_heat_map(total_rld, total_sample_data)

plotPCA(total_rld, intgroup=c("condition"))
plotPCA(total_rld, intgroup=c("sample"))

#####

results %<>% mutate(
  mean_fpkm = (AN1_fpkm + BN1_fpkm + CN1_fpkm +
                 AN2_fpkm + BN2_fpkm + CN2_fpkm) / 6)

expressed <- results %>% filter(mean_fpkm > 1)

ggplot(expressed, aes(l2fc_normal, l2fc_best)) +
  geom_point(colour="#444444") + 
  xlab(bquote(''~log[2]~'FC, normal STAR mapping')) + 
  ylab(bquote(''~log[2]~'FC, "best" filtering strategy')) + 
  theme_grey(base_size = 24)

cor(expressed$l2fc_normal, expressed$l2fc_best)

ggplot(expressed, aes(l2fc_normal, l2fc_conservative)) +
  geom_point(colour="#444444") + 
  xlab(bquote(''~log[2]~'FC, normal STAR mapping')) + 
  ylab(bquote(''~log[2]~'FC, "conservative" filtering strategy')) + 
  theme_grey(base_size = 24)

cor(expressed$l2fc_normal, expressed$l2fc_conservative)
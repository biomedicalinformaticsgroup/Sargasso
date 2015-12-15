library(DESeq2)
library(dplyr)
library(gplots)
library(magrittr)
library(RColorBrewer)
library(readr)
library(stringr)

# Read per-gene counts from a TSV file with two columns (Ensembl ID, read count)
read_counts <- function(count_file, sample_name, base_counts_dir="results/read_counts") {
  count_file_path = str_c(base_counts_dir, count_file, sep="/")
  counts <- read_tsv(count_file_path, col_names=c("gene", sample_name))
  return(counts)
}

# Set row names in a count data frame to the gene names, and remove the gene 
# name column
set_row_names <- function(count_data_frame) {
  row.names(count_data_frame) <- count_data_frame$gene
  count_data_frame %<>% select(-gene)
  return(count_data_frame)
}

# Perform differential expression with DESeq2
do_differential_expression <- function(count_data, sample_metadata, design_formula) {
  deseq_dataset <- DESeqDataSetFromMatrix(
    countData=count_data, colData=sample_metadata, design=design_formula)
  
  # Only retain genes with > 1 read across all samples
  deseq_dataset <- deseq_dataset[rowSums(counts(deseq_dataset)) > 1, ]
  
  total_dds <- DESeq(deseq_dataset)
}

# Plot a heatmap of sample-to-sample distances using regularized log 
# transformed count data
plot_sample_heatmap <- function(rld, sample_metadata) {
  total_distsRL <- dist(t(assay(rld)))
  total_mat <- as.matrix(total_distsRL)
  rownames(total_mat) <- colnames(total_mat) <- row.names(sample_metadata)
  total_hc <- hclust(total_distsRL)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(total_mat, Rowv=as.dendrogram(total_hc),
            symm=TRUE, trace="none",
            col = rev(hmcol), margin=c(13, 13))
}

# Read unfiltered and filtered per-gene counts - in this case from 3 mouse neuron 
# replicates (A, B, C) in two conditions (1, 2)
counts_AN1 <- read_counts("AN1.mouse.conventional.count", "AN1")
counts_AN1_filtered <- read_counts("AN1.mouse.count", "AN1_filtered")
counts_AN2 <- read_counts("AN2.mouse.conventional.count", "AN2")
counts_AN2_filtered <- read_counts("AN2.mouse.count", "AN2_filtered")
counts_BN1 <- read_counts("BN1.mouse.conventional.count", "BN1")
counts_BN1_filtered <- read_counts("BN1.mouse.count", "BN1_filtered")
counts_BN2 <- read_counts("BN2.mouse.conventional.count", "BN2")
counts_BN2_filtered <- read_counts("BN2.mouse.count", "BN2_filtered")
counts_CN1 <- read_counts("CN1.mouse.conventional.count", "CN1")
counts_CN1_filtered <- read_counts("CN1.mouse.count", "CN1_filtered")
counts_CN2 <- read_counts("CN2.mouse.conventional.count", "CN2")
counts_CN2_filtered <- read_counts("CN2.mouse.count", "CN2_filtered")

# Create a data frame containing unfiltered counts
count_data <- counts_AN1 %>%
  inner_join(counts_BN1) %>%
  inner_join(counts_CN1) %>%
  inner_join(counts_AN2) %>%
  inner_join(counts_BN2) %>%
  inner_join(counts_CN2)

# Create a data frame containing filtered counts
count_data_filtered <- counts_AN1_filtered %>%
  inner_join(counts_BN1_filtered) %>%
  inner_join(counts_CN1_filtered) %>%
  inner_join(counts_AN2_filtered) %>%
  inner_join(counts_BN2_filtered) %>%
  inner_join(counts_CN2_filtered)

# Create a data frame containing all counts, unfiltered and filtered
total_count_data <- count_data %>% inner_join(count_data_filtered)

# Set row names in count data frames to the gene names, as required for input to DESeq
count_data <- set_row_names(count_data)
count_data_filtered <- set_row_names(count_data_filtered)
total_count_data <- set_row_names(total_count_data)

# Create a data frame containing sample metadata for the unfiltered count data
sample_data <- data.frame(condition=c(rep("N1", 3), rep("N2", 3)),
                                sample=rep(c("A", "B", "C"), 2),
                                row.names=c("AN1", "BN1", "CN1",
                                            "AN2", "BN2", "CN2"))

# Create a data frame containing sample metadata for the filtered count data
sample_data_filtered <- data.frame(condition=c(rep("N1", 3), rep("N2", 3)),
                          sample=rep(c("A", "B", "C"), 2),
                          row.names=c("AN1_filtered", "BN1_filtered", "CN1_filtered",
                                      "AN2_filtered", "BN2_filtered", "CN2_filtered"))

# Create a data frame containing sample metadata for all the count data together
total_sample_data <-  data.frame(condition=c(rep("N1", 3), rep("N2", 3),
                                             rep("N1_filtered", 3), rep("N2_filtered", 3)),
                                 sample=rep(c("A", "B", "C"), 4),
                                 filtering=c(rep("No filtering", 6), rep("Filtering", 6)), 
                                 row.names=c("AN1", "BN1", "CN1",
                                             "AN2", "BN2", "CN2",
                                             "AN1_filtered", "BN1_filtered", "CN1_filtered",
                                             "AN2_filtered", "BN2_filtered", "CN2_filtered"))

# Create a DESeq dataset object for the total count data - this is only used 
# for diagnostic plots
total_dds <- do_differential_expression(
  total_count_data, total_sample_data, ~sample+condition)

# Calculate regularised log transformed counts
total_rld <- rlog(total_dds)

# Plot a heatmap of sample-to-sample distances for both unfiltered and filtered data
plot_sample_heatmap(total_rld, total_sample_data)

# PCA plots keys by experimental condition, sample number, and filtering state
plotPCA(total_rld, intgroup=c("condition"))
plotPCA(total_rld, intgroup=c("sample"))
plotPCA(total_rld, intgroup=c("filtering"))

# Create a DESeq dataset object for the unfiltered count data
dds <- do_differential_expression(count_data, sample_data, ~sample+condition)

# Get differential expression results for condition 2 vs 1 comparison 
# using unfiltered count data
res_2v1 <- results(dds, c("condition", "N2", "N1"))

# Create a DESeq dataset object for the unfiltered count data
dds_filtered <- do_differential_expression(
  count_data_filtered, sample_data_filtered, ~sample+condition)

# Get differential expression results for condition 2 vs 1 comparison 
# using filtered count data
res_2v1_filtered <- results(dds_filtered, c("condition", "N2", "N1"))

# Extract normalised counts from unfiltered data
counts <- dds %>% 
  counts(normalized=T) %>%
  as.data.frame() %>%
  add_rownames(var="gene")

# Extract raw un-normalised counts from unfiltered data
counts_raw <- dds %>% 
  counts(normalized=F) %>%
  as.data.frame() %>%
  add_rownames(var="gene") %>%
  rename(AN1_raw=AN1, BN1_raw=BN1, CN1_raw=CN1,
         AN2_raw=AN2, BN2_raw=BN2, CN2_raw=CN2)

# Extract normalised counts from filtered data
counts_filtered <- dds_filtered %>% 
  counts(normalized=T) %>%
  as.data.frame() %>%
  add_rownames(var="gene")

# Extract raw un-normalised counts from filtered data
counts_filtered_raw<- dds_filtered %>% 
  counts(normalized=F) %>%
  as.data.frame() %>%
  add_rownames(var="gene") %>%
  rename(AN1_filtered_raw=AN1_filtered, BN1_filtered_raw=BN1_filtered, CN1_filtered_raw=CN1_filtered,
         AN2_filtered_raw=AN2_filtered, BN2_filtered_raw=BN2_filtered, CN2_filtered_raw=CN2_filtered)

# Begin to construct a results object with DESeq results from both unfiltered 
# and filtered count data
results <- counts %>%
  inner_join(counts_raw) %>%
  inner_join(counts_filtered) %>%
  inner_join(counts_filtered_raw)

# Add columns for means of normalised counts in each condition for both 
# filtered and unfiltered data
results %<>%
  mutate(N1_mean_count=rowMeans(select(results, matches(".N1")))) %>%
  mutate(N2_mean_count=rowMeans(select(results, matches(".N2")))) %>%
  mutate(N1_mean_count_filtered=rowMeans(select(results, matches(".N1_filtered")))) %>%
  mutate(N2_mean_count_filtered=rowMeans(select(results, matches(".N2_filtered"))))

# Add column containing the per-gene fraction of reads lost across all samples due to
# filtering
results %<>%
  mutate(mean_frac_lost=1 - (AN1_filtered_raw + BN1_filtered_raw + CN1_filtered_raw + 
                             AN2_filtered_raw + BN2_filtered_raw + CN2_filtered_raw) / 
                               (AN1_raw + BN1_raw + CN1_raw + AN2_raw + BN2_raw + CN2_raw))

# Add differential expression results columns for unfiltered data:
# - raw log2 fold change calculated from mean normalised counts
# - DESeq2 moderated fold change
# - DESeq2 p-value
# - DESeq2 adjusted p-value
results %<>% mutate(raw_l2fc_2v1=log2(N2_mean_count/N1_mean_count))
results %<>% inner_join(res_2v1 %>% as.data.frame() %>% add_rownames(var="gene"), by="gene") %>%
  rename(l2fc_2v1 = log2FoldChange, pval_2v1 = pvalue, padj_2v1 = padj) %>%
  select(-baseMean, -lfcSE, -stat)

# Add the same differential expression results columns for the filtered data
results %<>% mutate(raw_l2fc_2v1_filtered=log2(N2_mean_count_filtered/N1_mean_count_filtered))
results %<>% inner_join(res_2v1_filtered %>% as.data.frame() %>% add_rownames(var="gene"), by="gene") %>%
  rename(l2fc_2v1_filtered = log2FoldChange, pval_2v1_filtered = pvalue, padj_2v1_filtered = padj) %>%
  select(-baseMean, -lfcSE, -stat)

# Add columns containing gene ranks according to DESeq log2 fold change and p-value in
# both the filtered and unfiltered data
results %<>% arrange(desc(l2fc_2v1)) %>% mutate(l2fc_2v1_rank=row_number())
results %<>% arrange(desc(l2fc_2v1_filtered)) %>% mutate(l2fc_2v1_rank_filtered=row_number())
results %<>% arrange(pval_2v1) %>% mutate(pval_2v1_rank=row_number())
results %<>% arrange(pval_2v1_filtered) %>% mutate(pval_2v1_rank_filtered=row_number())

# Read in information about mouse genes, downloaded from Ensembl BioMart via the
# following query:
#
#<?xml version="1.0" encoding="UTF-8"?>
#  <!DOCTYPE Query>
#  <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" > 
#  <Dataset name = "mmusculus_gene_ensembl" interface = "default" >
#  <Attribute name = "ensembl_gene_id" />
#  <Attribute name = "chromosome_name" />
#  <Attribute name = "external_gene_name" />
#  <Attribute name = "description" />
#  </Dataset>
#  </Query>
mouse_genes <- read_tsv("data/ensembl-82/mouse_genes.tsv", 
                        col_types = list(`Chromosome Name` = col_character())) %>%
  dplyr::rename("gene"=`Ensembl Gene ID`,
                "gene_name"=`Associated Gene Name`,
                "chromosome"=`Chromosome Name`,
                "description"=Description) 

# Add gene info to results object
results %<>% left_join(mouse_genes)

# Order results data table for output...
results %<>% select(gene, gene_name, description, chromosome,
                    AN1, BN1, CN1, N1_mean_count, 
                    AN2, BN2, CN2, N2_mean_count,
                    AN1_raw, BN1_raw, CN1_raw, 
                    AN2_raw, BN2_raw, CN2_raw,
                    AN1_filtered, BN1_filtered, CN1_filtered, N1_mean_count_filtered,
                    AN2_filtered, BN2_filtered, CN2_filtered, N2_mean_count_filtered,
                    AN1_filtered_raw, BN1_filtered_raw, CN1_filtered_raw, 
                    AN2_filtered_raw, BN2_filtered_raw, CN2_filtered_raw,
                    raw_l2fc_2v1, l2fc_2v1, l2fc_2v1_rank, pval_2v1, pval_2v1_rank, padj_2v1,
                    raw_l2fc_2v1_filtered, l2fc_2v1_filtered, l2fc_2v1_rank_filtered, pval_2v1_filtered, pval_2v1_rank_filtered, padj_2v1_filtered,
                    mean_frac_lost)

# ...and write to a CSV file
write_csv(results, "results/differential_expression/diff_expr.csv")

# Filter results to retain only those with mean normalised count > 100. Add ranks
# by DESeq p-value in both filtered and unfiltered data. ("expressed" genes)
results_expr <- results %>% filter((N1_mean_count + N2_mean_count) / 2 > 100)
results_expr %<>% arrange(pval_2v1) %>% mutate(pval_2v1_rank_restricted=row_number())
results_expr %<>% arrange(pval_2v1_filtered) %>% mutate(pval_2v1_rank_filtered_restricted=row_number())

# Filter results to retain only those with adjusted p-value < 0.1 (in *unfiltered* 
# DESeq results). Add ranks byDESeq p-value in both filtered and unfiltered data. 
# ("significant" genes)
sig_results <- results %>% filter(padj_2v1 < 0.1)
sig_results %<>% arrange(pval_2v1) %>% mutate(pval_2v1_rank_restricted=row_number())
sig_results %<>% arrange(pval_2v1_filtered) %>% mutate(pval_2v1_rank_filtered_restricted=row_number())

# Filter results to retain only those with adjusted p-value < 0.1 *and* mean 
# normalised count > 100. Add ranks by DESeq p-value in both filtered and
# unfiltered data. ("expressed and significant" genes)
sig_results_expr <- sig_results %>% filter((N1_mean_count + N2_mean_count) / 2 > 100)
sig_results_expr %<>% arrange(pval_2v1) %>% mutate(pval_2v1_rank_restricted=row_number())
sig_results_expr %<>% arrange(pval_2v1_filtered) %>% mutate(pval_2v1_rank_filtered_restricted=row_number())

# Filter results to retain only those with adjusted p-value < 0.1 in *filtered* data 
sig_results_filtered <- results %>% filter(padj_2v1_filtered < 0.1)

# Further filter results to retain only those with adjusted p-value < 0.1 in *filtered*
# data and mean normalised count > 100 in *unfiltered* data
sig_results_expr_filtered <- sig_results_filtered %>% filter((N1_mean_count + N2_mean_count) / 2 > 100)

# Set some graphics parameters for following scatter plots
par(pch=20, col=rgb(0, 0, 0, 0.3))

# Scatter plot of log2 fold change in filtered vs unfiltered DESeq results, all genes

plot(results$l2fc_2v1, results$l2fc_2v1_filtered, 
     xlab="Log2 fold change", 
     ylab="Log2 fold change (filtered)",
     xlim=c(-4, 9), ylim=c(-4, 9))

# Scatter plot of log2 fold change in filtered vs unfiltered DESeq results, 
# "significant" genes only

plot(sig_results$l2fc_2v1, sig_results$l2fc_2v1_filtered, 
     xlab="Log2 fold change", 
     ylab="Log2 fold change (filtered)",
     xlim=c(-4, 9), ylim=c(-4, 9))

# Scatter plot of log2 fold change in filtered vs unfiltered DESeq results, 
# "expressed" genes only

plot(results_expr$l2fc_2v1, results_expr$l2fc_2v1_filtered, 
     xlab="Log2 fold change", 
     ylab="Log2 fold change (filtered)",
     xlim=c(-4, 9), ylim=c(-4, 9))

# Scatter plot of log2 fold change in filtered vs unfiltered DESeq results,
# "significant expressed" genes only

plot(sig_results_expr$l2fc_2v1, sig_results_expr$l2fc_2v1_filtered, 
     xlab="Log2 fold change", 
     ylab="Log2 fold change (filtered)",
     xlim=c(-4, 9), ylim=c(-4, 9))

# Scatter plot of log10 p-value in filtered vs unfiltered results, all genes

plot(log10(results$pval_2v1), log10(results$pval_2v1_filtered),
     xlab="Log10 raw p-value", ylab="Log10 raw p-value (filtered)",
     xlim=c(-300, 0), ylim=c(-300, 0))

# Scatter plot of log10 p-value in filtered vs unfiltered results, "significant"
# genes only

plot(log10(sig_results$pval_2v1), log10(sig_results$pval_2v1_filtered),
     xlab="Log10 raw p-value", ylab="Log10 raw p-value (filtered)",
     xlim=c(-300, 0), ylim=c(-300, 0))

# Scatter plot of log10 p-value in filtered vs unfiltered results, "expressed"
# genes only

plot(log10(results_expr$pval_2v1), log10(results_expr$pval_2v1_filtered),
     xlab="Log10 raw p-value", ylab="Log10 raw p-value (filtered)",
     xlim=c(-300, 0), ylim=c(-300, 0))

# Scatter plot of log10 p=value in filtered vs unfiltered results, 
# "significant expressed" genes only

plot(log10(sig_results_expr$pval_2v1), log10(sig_results_expr$pval_2v1_filtered),
     xlab="Log10 raw p-value", ylab="Log10 raw p-value (filtered)",
     xlim=c(-300, 0), ylim=c(-300, 0))

# Scatter plot of p-value rank in filtered vs unfiltered results, all genes. Vertical
# line indicates FDR < 0.1 cutoff.

plot(results$pval_2v1_rank, results$pval_2v1_rank_filtered,
     xlab="p-value rank", ylab="p-value rank (filtered)",
     xlim=c(0, 25500), ylim=c(0, 25500))
abline(v=nrow(sig_results))

# Scatter plot of p-value rank in filtered vs unfiltered results, "significant" genes,
# rank amongst *all* genes

plot(sig_results$pval_2v1_rank, sig_results$pval_2v1_rank_filtered,
     xlab="p-value rank", ylab="p-value rank (filtered)",
     ylim=c(0, 25000))

# Scatter plot of p-value rank in filtered vs unfiltered results, "significant" genes,
# rank amongst "signficant" genes

plot(sig_results$pval_2v1_rank_restricted, sig_results$pval_2v1_rank_filtered_restricted,
     xlab="p-value rank", ylab="p-value rank (filtered)")

# Scatter plot of p-value rank in filtered vs unfiltered results, "expressed" genes,
# rank amongst all genes. Vertical line indicates FDR < 0.1 cutoff.

plot(results_expr$pval_2v1_rank, results_expr$pval_2v1_rank_filtered,
     xlab="p-value rank", ylab="p-value rank (filtered)",
     xlim=c(0, 25500), ylim=c(0, 25500))
abline(v=nrow(sig_results))

# Scatter plot of p-value rank in filtered vs unfiltered results, "expressed" genes,
# rank amongst "expressed" genes. Vertical line indicates FDR < 0.1 cutoff.

plot(results_expr$pval_2v1_rank_restricted, results_expr$pval_2v1_rank_filtered_restricted,
     xlab="p-value rank", ylab="p-value rank (filtered)")
abline(v=nrow(sig_results_expr))

# Scatter plot of p-value rank in filtered vs unfiltered results, 
# "significant expressed" genes, rank amongst *all* genes

plot(sig_results_expr$pval_2v1_rank, sig_results_expr$pval_2v1_rank_filtered,
     xlab="p-value rank", ylab="p-value rank (filtered)",
     ylim=c(0, 25000))

# Scatter plot of p-value rank in filtered vs unfiltered results, 
# "significant expressed" genes, rank amongst "signficant expressed" genes

plot(sig_results_expr$pval_2v1_rank_restricted, sig_results_expr$pval_2v1_rank_filtered_restricted,
     xlab="p-value rank", ylab="p-value rank (filtered)")

# Plot of difference in p-value rank in filtered vs unfiltered results against 
# mean normalised count, all genes

plot(log10(results$N1_mean_count + results$N2_mean_count)/2, results$pval_2v1_rank - results$pval_2v1_rank_filtered,
     xlab="log10 mean normalised count", ylab="p-value rank difference")

# Plot of difference in p-value rank in filtered vs unfiltered results against 
# mean normalised count, "significant" genes

plot(log10(sig_results$N1_mean_count + sig_results$N2_mean_count)/2, sig_results$pval_2v1_rank - sig_results$pval_2v1_rank_filtered,
     xlab="log10 mean normalised count", ylab="p-value rank difference")

# Plot of difference in log2 fold change in filtered vs unfiltered results
# against mean normalised count, all genes

plot(log10(results$N1_mean_count + results$N2_mean_count)/2, results$l2fc_2v1 - results$l2fc_2v1_filtered,
     xlab="log10 mean normalised count", ylab="log2 FC difference")

# Plot of difference in log2 fold change in filtered vs unfiltered results,
# against mean normalised count, "significant" genes

plot(log10(sig_results$N1_mean_count + sig_results$N2_mean_count)/2, sig_results$l2fc_2v1 - sig_results$l2fc_2v1_filtered,
     xlab="log10 mean normalised count", ylab="log2 FC difference")

# Plot of difference in p-value rank in filtered vs unfiltered results against
# fraction of reads lost due to filtering, all genes

plot(results$mean_frac_lost, results$pval_2v1_rank - results$pval_2v1_rank_filtered,
     xlab="Mean fraction of reads lost in filtering", ylab="p-value rank difference",
     xlim=c(0, 1), ylim=c(-17000, 17000))

# Plot of difference in p-value rank in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "significant" genes

plot(sig_results$mean_frac_lost, sig_results$pval_2v1_rank - sig_results$pval_2v1_rank_filtered,
     xlab="Mean fraction of reads lost in filtering", ylab="p-value rank difference",
     xlim=c(0, 1), ylim=c(-17000, 17000))

# Plot of difference in p-value rank in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "significant" genes (rank amongst 
# "significant" genes)

plot(sig_results$mean_frac_lost, sig_results$pval_2v1_rank_restricted - sig_results$pval_2v1_rank_filtered_restricted,
     xlab="Mean fraction of reads lost in filtering", ylab="p-value rank difference",
     xlim=c(0, 1))

# Plot of difference in p-value rank in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "expressed" genes

plot(results_expr$mean_frac_lost, results_expr$pval_2v1_rank - results_expr$pval_2v1_rank_filtered,
     xlab="Mean fraction of reads lost in filtering", ylab="p-value rank difference",
     xlim=c(0, 1), ylim=c(-17000, 17000))

# Plot of difference in p-value rank in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "expressed" genes (rank amongst 
# "expressed" genes)

plot(results_expr$mean_frac_lost, results_expr$pval_2v1_rank_restricted - results_expr$pval_2v1_rank_filtered_restricted,
     xlab="Mean fraction of reads lost in filtering", ylab="p-value rank difference",
     xlim=c(0, 1))

# Plot of difference in p-value rank in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "significant expressed" genes

plot(sig_results_expr$mean_frac_lost, sig_results_expr$pval_2v1_rank - sig_results_expr$pval_2v1_rank_filtered,
     xlab="Mean fraction of reads lost in filtering", ylab="p-value rank difference",
     xlim=c(0, 1), ylim=c(-17000, 17000))

# Plot of difference in p-value rank in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "significant expressed" genes (rank amongst 
# "significant expressed" genes)

plot(sig_results_expr$mean_frac_lost, sig_results_expr$pval_2v1_rank_restricted - sig_results_expr$pval_2v1_rank_filtered_restricted,
     xlab="Mean fraction of reads lost in filtering", ylab="p-value rank difference",
     xlim=c(0, 1))

# Plot of difference in log2 fold change in filtered vs unfiltered results against
# fraction of reads lost due to filtering, all genes

plot(results$mean_frac_lost, results$l2fc_2v1 - results$l2fc_2v1_filtered,
     xlab="Mean fraction of reads lost in filtering", ylab="log2 fold change difference",
     xlim=c(0, 1), ylim=c(-3.2, 3.2))

# Plot of difference in log2 fold change in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "significant" genes

plot(sig_results$mean_frac_lost, sig_results$l2fc_2v1 - sig_results$l2fc_2v1_filtered,
     xlab="Mean fraction of reads lost in filtering", ylab="log2 fold change difference",
     xlim=c(0, 1), ylim=c(-3.2, 3.2))

# Plot of difference in log2 fold change in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "expressed" genes

plot(results_expr$mean_frac_lost, results_expr$l2fc_2v1 - results_expr$l2fc_2v1_filtered,
     xlab="Mean fraction of reads lost in filtering", ylab="log2 fold change difference",
     xlim=c(0, 1), ylim=c(-3.2, 3.2))

# Plot of difference in log2 fold change in filtered vs unfiltered results against
# fraction of reads lost due to filtering, "significant expressed" genes

plot(sig_results_expr$mean_frac_lost, sig_results_expr$l2fc_2v1 - sig_results_expr$l2fc_2v1_filtered,
     xlab="Mean fraction of reads lost in filtering", ylab="log2 fold change difference",
     xlim=c(0, 1), ylim=c(-3.2, 3.2))

# Total number of genes in DESeq analysis
nrow(results)

# Number of "expressed" genes
nrow(results_expr)

# Number of "significant" genes
nrow(sig_results)

# Number of "significant" genes for which adjusted p-value in filtered results is also < 0.1
nrow(sig_results %>% filter(padj_2v1_filtered < 0.1))

# Number of "significant expressed" genes
nrow(sig_results_expr)

# Number of "significant expressed" genes for which adjusted p-value in filtered results
# is also < 0.1
nrow(sig_results_expr %>% filter(padj_2v1_filtered < 0.1))

# Number of genes with FDR < 0.1 in *filtered* data
nrow(sig_results_filtered)

# Number of genes with FDR < 0.1 in *filtered* data for which mean normalised 
# count > 100 in *unfiltered* data
nrow(sig_results_expr_filtered)

# Correlations between DESeq results in filtered and unfiltered data

# Correlation in p-value rank, all genes
cor(results$pval_2v1_rank, results$pval_2v1_rank_filtered)

# Correlation in log2 fold change, all genes
cor(results$l2fc_2v1, results$l2fc_2v1_filtered)

# Correlation in p-value rank, "expressed" genes
cor(results_expr$pval_2v1_rank, results_expr$pval_2v1_rank_filtered)

# Correlation in log2 fold change, "expressed" genes
cor(results_expr$l2fc_2v1, results_expr$l2fc_2v1_filtered)

# Correlation in p-value rank, "significant" genes
cor(sig_results$pval_2v1_rank, sig_results$pval_2v1_rank_filtered)

# Correlation in log2 fold change, "significant" genes
cor(sig_results$l2fc_2v1, sig_results$l2fc_2v1_filtered)

# Correlation in p-value rank, "significant expressed" genes
cor(sig_results_expr$pval_2v1_rank, sig_results_expr$pval_2v1_rank_filtered)

# Correlation in log2 fold change "significant expressed" genes
cor(sig_results_expr$l2fc_2v1, sig_results_expr$l2fc_2v1_filtered)

##########

#results_mfl <- results %>% filter(mean_frac_lost < 0.5)
#results_expr_mfl <- results_expr %>% filter(mean_frac_lost < 0.5)

#sig_results_mfl <- sig_results %>% filter(mean_frac_lost < 0.5)
#sig_results_expr_mfl <- sig_results_expr %>% filter(mean_frac_lost < 0.5)

#sig_results_mfl_filtered <- sig_results_filtered %>% filter(mean_frac_lost < 0.5)
#sig_results_expr_mfl_filtered <- sig_results_expr_filtered %>% filter(mean_frac_lost < 0.5)

#nrow(results_mfl)
#nrow(results_expr_mfl)

#cor(results_mfl$pval_2v1_rank, results_mfl$pval_2v1_rank_filtered)
#cor(results_mfl$l2fc_2v1, results_mfl$l2fc_2v1_filtered)
#cor(results_expr_mfl$pval_2v1_rank, results_expr_mfl$pval_2v1_rank_filtered)
#cor(results_expr_mfl$l2fc_2v1, results_expr_mfl$l2fc_2v1_filtered)

#cor(sig_results_mfl$pval_2v1_rank, sig_results_mfl$pval_2v1_rank_filtered)
#cor(sig_results_mfl$l2fc_2v1, sig_results_mfl$l2fc_2v1_filtered)
#cor(sig_results_expr_mfl$pval_2v1_rank, sig_results_expr_mfl$pval_2v1_rank_filtered)
#cor(sig_results_expr_mfl$l2fc_2v1, sig_results_expr_mfl$l2fc_2v1_filtered)

#nrow(sig_results_mfl)
#nrow(sig_results_mfl %>% filter(padj_2v1_filtered < 0.1))
#nrow(sig_results_expr_mfl)
#nrow(sig_results_expr_mfl %>% filter(padj_2v1_filtered < 0.1))

#nrow(sig_results_mfl_filtered)
#nrow(sig_results_mfl_filtered %>% filter(padj_2v1 < 0.1))
#nrow(sig_results_expr_mfl_filtered)
#nrow(sig_results_expr_mfl_filtered %>% filter(padj_2v1 < 0.1))

#plot(results_expr_mfl$l2fc_2v1, results_expr_mfl$l2fc_2v1_filtered, 
#     xlab="Log2 fold change", 
#     ylab="Log2 fold change (filtered)",
#     xlim=c(-4, 9), ylim=c(-4, 9))

#plot(log10(results_expr_mfl$pval_2v1), log10(results_expr_mfl$pval_2v1_filtered),
#     xlab="Log10 raw p-value", ylab="Log10 raw p-value (filtered)",
#     xlim=c(-300, 0), ylim=c(-300, 0))

#plot(results_expr_mfl$pval_2v1_rank, results_expr_mfl$pval_2v1_rank_filtered,
#     xlab="p-value rank", ylab="p-value rank (filtered)",
#     xlim=c(0, 25500), ylim=c(0, 25500))
#abline(v=nrow(sig_results))

#plot(log10(results$N1_mean_count_filtered + results$N2_mean_count_filtered)/2, results$pval_2v1_rank - results$pval_2v1_rank_filtered,
#     xlab="log10 mean normalised count (filtered)", ylab="p-value rank difference")

#plot(log10(sig_results$N1_mean_count_filtered + sig_results$N2_mean_count_filtered)/2, sig_results$pval_2v1_rank - sig_results$pval_2v1_rank_filtered,
#     xlab="log10 mean normalised count (filtered)", ylab="p-value rank difference")

#plot(log10(results$N1_mean_count_filtered + results$N2_mean_count_filtered)/2, results$l2fc_2v1 - results$l2fc_2v1_filtered,
#     xlab="log10 mean normalised count (filtered)", ylab="log2 FC difference")

#plot(log10(sig_results$N1_mean_count_filtered + sig_results$N2_mean_count_filtered)/2, sig_results$l2fc_2v1 - sig_results$l2fc_2v1_filtered,
#     xlab="log10 mean normalised count (filtered)", ylab="log2 FC difference")

#plot(log10(results_expr_mfl$N1_mean_count + results_expr_mfl$N2_mean_count)/2, results_expr_mfl$pval_2v1_rank - results_expr_mfl$pval_2v1_rank_filtered,
#     xlab="log10 mean normalised count", ylab="p-value rank difference")

#plot(log10(results_expr_mfl$N1_mean_count + results_expr_mfl$N2_mean_count)/2, results_expr_mfl$l2fc_2v1 - results_expr_mfl$l2fc_2v1_filtered,
#     xlab="log10 mean normalised count", ylab="log2 FC difference")

#plot(results_expr_mfl$mean_frac_lost, results_expr_mfl$pval_2v1_rank - results_expr_mfl$pval_2v1_rank_filtered,
#     xlab="Mean fraction of reads lost in filtering", ylab="p-value rank difference",
#     xlim=c(0, 1), ylim=c(-17000, 17000))

#plot(results_expr_mfl$mean_frac_lost, results_expr_mfl$l2fc_2v1 - results_expr_mfl$l2fc_2v1_filtered,
#     xlab="Mean fraction of reads lost in filtering", ylab="log2 fold change difference",
#     xlim=c(0, 1), ylim=c(-3.2, 3.2))

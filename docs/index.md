Welcome to Sargasso
===================

*Sargasso* is a tool to separate mixed-species RNA-seq reads according to their species of origin.

[RNA-sequencing](http://en.wikipedia.org/wiki/RNA-Seq) has become an important technique in cellular biology for characterising and quantifying the [transcriptomes](http://en.wikipedia.org/wiki/Transcriptome) of particular species, and for analysing the [differential expression](https://en.wikipedia.org/wiki/RNA-Seq#Differential_expression_and_absolute_quantification_of_transcripts) of genes and transcripts. However, a number of recent experimental techniques produce RNA-seq data originating from a mixture of species. Previously, researchers have developed ad-hoc solutions to separate such RNA-seq reads into species-specific sets; *Sargasso* is an efficient, reliable tool to perform this task, achieving high specificity and sensitivity, while requiring minimal setup and intervention by the user.

Given an input set of [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format>) files for a number of samples, each of which contain mixed-species RNA-seq read data, *Sargasso* separates reads according to their true species of origin, and outputs per-sample [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files describing the mapping of these reads to the respective genomes. While the tool allows the user fine control at each stage, in normal usage the whole pipeline can be executed automatically with a single command.

For further details, please see below:

1. [Installation](installation.md)
2. [Example usage](example_usage.md)
3. [Pipeline description](pipeline.md)
4. [Usage reference](usage_reference.md)
5. [Support scripts](support_scripts.md)
6. [References](references.md)

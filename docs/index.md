*Sargasso*<sup>*</sup> is a tool to separate mixed-species high-throughput sequencing reads according to their species of origin.

High-throughput sequencing techniques, such as [RNA-seq](https://en.wikipedia.org/wiki/RNA-Seq), [ChIP-seq](https://en.wikipedia.org/wiki/ChIP-sequencing) and [ATAC-seq](https://en.wikipedia.org/wiki/ATAC-seq), have become important techniques in cellular biology, for characterising and quantifying the transcriptome of particular species, for analysing protein interactions with their DNA, and for assessing chromatin accessibility. However, a number of recent experimental techniques produce high-throughput sequencing data originating from a mixture of species. Often, researchers have developed ad-hoc solutions to separate such sequencing reads into species-specific sets; *Sargasso* is an efficient, reliable tool to perform this task, achieving high specificity and sensitivity even for closely-related species, while requiring minimal setup and intervention by the user.

Given an input set of [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files for a number of samples, each of which contain mixed-species sequencing read data, *Sargasso* separates reads according to their true species of origin, and outputs per-sample [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files describing the mapping of these reads to the respective genomes. While the tool allows the user fine control at each stage, in normal usage the whole pipeline can be executed automatically with a single command.

For further details, please see below:

1. [Installation](installation.md)
2. [Example usage](example_usage.md)
3. [Pipeline description](pipeline.md)
4. [Usage reference](usage_reference.md)
5. [Support scripts](support_scripts.md)
6. [References](references.md)

<sup>* **S**argasso **A**ssigns **R**eads to **G**enomes **A**ccording to **S**pecies-**S**pecific **O**rigin</sup>

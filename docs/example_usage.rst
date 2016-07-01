Example usage
=============

To illustrate the usage of the *species-separator* pipeline, we will process a test data set consisting of RNA-seq reads derived from both mouse and rat. The test FASTQ files can be found in the directory ``pipeline_test/data/fastq/`` within the *species-separator* repository.

We begin by constructing a tab-separated file giving the paths of the raw RNA-seq read data files for each sample in our experiment. This should contain one line per-sample, giving a sample name, and two comma-separated lists of FASTQ files containing paired-end RNA-seq reads (or a single comma-separated list in the case of single-end reads). Note that the FASTQ files are assumed to be gzipped.

In our example case we have a single sample, with a single pair of paired-end read files, and hence our ``test_samples.tsv`` file is particularly simple. It contains the single line::

    our_sample  pipeline_test/data/fastq/mouse_rat_test_1.fastq.gz  pipeline_test/data/fastq/mouse_rat_test_2.fastq.gz

The *species-separator* pipeline uses STAR [STAR]_, an efficient and accurate short RNA-seq read aligner, to map reads to reference genomes. We will assume that STAR indices have already been built for the mouse and rat genomes, and are located in the directories ``~/data/genome/<species>/STAR_index/``. Then the entire species separation pipeline can be executed with the following command::

    species_separator --reads-base-dir=<repository_dir> --s1-index=~/data/genome/mouse/STAR_index --s2-index=~/data/genome/rat/STAR_index --best --run-separation mouse rat test_samples.tsv test_results

where ``<repository_dir>`` is the full path to the repository directory.

This command will execute the species separation pipeline in the background, using ``nohup`` to run commands immune from hangup. Results are output to directory ``test_results``, and pipeline progress can be monitored by examining the file ``nohup.out`` that is written in this directory.

On this small data set (100,000 paired-end reads), species separation should take a matter of minutes. On completion, the ``test_results`` directory will contain the following sub-directories:

* ``star_indices``: Contains a link to the STAR index directory for each species.
* ``raw_reads``: Contains links to the gzipped FASTQ files for each sample.
* ``mapped_reads``: Contains a BAM file for each sample and species, describing the mapping of the RNA-seq reads in that sample to the species' genome.
* ``sorted_reads``: Contains a BAM file for each sample and species, where the mapped reads above have been sorted in name order.
* ``filtered_reads``: Contains a BAM file for each sample and species, describing the mapping of the RNA-seq reads determined to have originated from that species, to that species' genome.

The BAM files in the ``filtered_reads`` directory are the final output of the *species-separator* pipeline. These can then be taken as input to further downstream analyses, for example for read counting and differential expression.

In addition, two further log files are written. In the ``filtered_reads`` directory, ``overall_filtering_summary.txt`` contains per-sample statistics describing the reads that were assigned to each genome, or were rejected as ambiguous. In the top-level ``test_results`` directory, ``execution_record.txt`` contains a record of the command line options that were passed to *species-separator*, and the date and time of execution.

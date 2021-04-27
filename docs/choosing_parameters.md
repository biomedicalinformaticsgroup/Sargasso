Choosing parameters
===================

While Sargasso includes a number of "pre-packaged" filtering strategies, each providing a particular balance between sensitivity and specificity, users may be interested to explore choices of filtering parameter in greater depth, allowing them to choose parameters that are tailored to their own experimental situation and the sequencing technology used. The ``sargasso_parameter_test`` script is provided to aid this process.

Usage
=====

The ``sargasso_parameter_test`` script runs Sargasso a number of times, on a set of provided samples, using a different filtering parameters each time. Filtering statistics are gathered, and informative plots of sensitivity and specificity are produced, allowing users to examine the trade-offs between these two measures. 

n.b. it is intended that ``sargasso_parameter_test`` be run on **single-species samples**, thus allowing incorrectly assigned reads to be easily identified.

    sargasso_parameter_test.sh <data-type>
        [--mapper-executable <mapper-executable>]
        [--reads-base-dir <reads-base-dir>]
        [--mismatch-setting '<mismatch-values>']
        [--minmatch-setting '<minmatch-values>']
        [--multimap-setting '<multimap-values>']
        [--num-threads <num-threads>]
        [--plot-format <plot-format>]
        [--skip-init-run]
        --samples-origin '<species> <species> ...'
        <samples-file> <output-dir>
        (<species> <species-info>)
        (<species> <species-info>)
         ...

* ``<data-type>`` (_text parameter_): One of 'dnaseq' or 'rnaseq' (matches Sargasso parameter of same name).
* ``--mapper-executable <mapper-executable>`` (_file path_): Alignment tool executable path (matches Sargasso parameter).
* ``--reads-base-dir <reads-base-dir>`` (_file path_): Base directory for raw sequencing files (matches Sargasso parameter).
* ``--mismatch-setting '<mismatch-values>'`` (_list of integers_): List of values for maximum percentage of read bases allowed to be mismatched against the genome during filtering.
* ``--minmatch-setting '<minmatch-values>'`` (_list of integers_): List of values for maximum percentage of read length allowed to not be mapped during filtering.
* ``--multimap-setting '<multimap-values>'`` (_list of integers_): List of values for maximum number of multiple mappings allowed during filtering.
* ``--num-threads <num-threads>`` (_integer_): Number of threads for parallel processing (matches Sargasso parameter).
* ``--plot-format <plot-format>`` (_text parameter_): One of 'pdf' or 'png'.
* ``--skip-init-run`` (_flag_): If set, do not run initial mapping step (it is assumed that this has already been run).
* ``--samples-origin '<species> <species> ...'`` (_list of text parameters_): List of true species of origin of the input samples (one for each sample).
* ``<samples-file>`` (_file path_): TSV file giving paths of raw sequencing read data files for each sample (matches Sargasso parameter).
* ``<output-dir>`` (_file path_): Output directory into which results will be written.
* ``<species>`` (_text parameter_): Name of the nth species that the input samples are to be separated into (matches Sargasso parameter).
* ``<species-info>`` (_text parameter_): Alignment tool for the nth speciees that the input samples are to be separated into (matches Sargasso parameter).

Output
======

Sargasso filtering will be run for every combination of the values of the parameters ``--mismatch-setting``, ``--minmatch-setting`` and ``--multimap-setting``, though initial mapping will only be run once (or not at all, if the flag ``--skip-init-run`` is set, in which case mapping is assumed to have already been run using this script). The main output of ``sargasso_parameter_test`` are two plots – ``counts.{png|pdf}`` and ``percentage.{png|pdf}``. These show, respectively, the number or percentage of reads which are assigned, marked as ambiguous, or rejected, for each (i) sample, (ii) species into which data is being separated, and (iii) combination of the mismatch, minmatch and multimap parameters.

Examining these graphs can aid the user to choose what filtering parameter values can be used to achieve an acceptable assignment of data to the correct species of origin, while minimising misassignment to the wrong species. 
Example usage
=============

    sargasso_parameter_test rnaseq 
        --samples-origin 'rat' 
        --mismatch-setting '0 2 4' 
        --minmatch-setting '0 2 4' 
        --multimap-setting '1' 
        --plot-format png 
        test_sample.tsv ~/tmp/sargasso_test 
        rat /srv/data/genome/rat/ensembl-103/STAR_indices/toplevel 
        human /srv/data/genome/human/ensembl-103/STAR_indices/primary_assembly

In this case, we are running the test script using a single sample, whose true species of origin is rat (paths to the raw reads files for the sample are contained in the file ``test_sample.tsv``). Sargasso filtering will be run for every combination of the specified mismatch (0, 2, 4), minmatch (0, 2, 4), and multimap (1) settings – a total of 9 runs – separating the input sample in rat and human reads (where any reads assigned to human will be incorrect assignments). The output files ``counts.png`` and ``percentage.png`` will be written into the ouptu directory ``~/tmp/sargasso_test``.

Usage reference
===============

The main ``species_separator`` script has a number of command line options to alter the behaviour and function of read assignment:

    species_separator <data-type>
        [--log-level=<log-level>]
        [--reads-base-dir=<reads-base-dir>] [--num-threads=<num-threads>]
        [--mismatch-threshold=<mismatch-threshold>]
        [--minmatch-threshold=<minmatch-threshold>]
        [--multimap-threshold=<multimap-threshold>]
        [--reject-multimaps] 
        [--best] [--conservative] [--recall] [--permissive]
        [--run-separation]
        [--delete-intermediate]
        [--mapper-executable=<mapper-executable>]
        [--mapper-index-executable=<mapper-index-executable>]
        [--sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>]
        <samples-file> <output-dir>
        (<species> <species-info>)
        (<species> <species-info>)
        ...

These parameters are listed and explained below for reference, grouped by their functionality.

Core
----

These parameters are required as the base minimum for the execution of the pipeline.

* ``<data-type>`` (_text parameter_): One of "dnaseq" or "rnaseq".
* ``<samples-file>`` (_file path_): TSV file giving paths (relative to ``<reads-base-dir>``) of raw sequencing read data files for each sample.
* ``<output-dir>`` (_file path_): Output directory into which the Makefile will be written, and in which species separation will be performed.
* ``<species>`` (_text parameter_): Name of the nth species. Note that at least two species must be specified. While there is no theoretical upper bound on the number of species, the depth of sequencing required to obtain enough reads for each species after separation may impose a practical limit.
* ``<species-info>`` (_file paths_): Alignment tool information for the nth species. For DNA sequencing data, this parameter should either consist of (i) the path to a Bowtie2 index directory for the species _or_ (ii) a FASTA file containing genome sequences for the species. For RNA sequencing data, the parameter should either consist of (i) the path to the STAR index directory for the species _or_ (ii) a comma separated list of two elements, the first of which is the path to a GTF annotation file for the species, and the second is the path to a directory containing genome FASTA files for the species.
    
Mapping
-------

These parameters control the mapping of mixed-species RNA-seq reads to reference genomes.

* ``--reads-base-dir=<reads-base-dir>`` (_file path_): Base directory for raw RNA-seq read data files.
* ``--mapper-executable`` (_file path_): Specifies the alignment tool executable path --- use this option to run Sargasso with a particular version of either Bowtie2 or STAR.
* ``--mapper-index-executable`` (_file path_): For DNA sequencing data, specifies the Bowtie2 index building tool path --- use this option to run Sargasso with a particular version of ``bowtie2-build`` (n.b. for RNA sequencing data, this option is ignored).

Assignment criteria and optimisation
----------------------------------

These parameters are used to specify criteria that affect how reads are assigned to each species, either by choosing individual values for each parameter or by selecting pre-configured assignment profiles.

* ``--mismatch-threshold=<mismatch-threshold>`` (_float_): Maximum percentage of read bases allowed to be mismatches against the genome during filtering (default: 0).
* ``--minmatch-threshold=<minmatch-threshold>`` (_float_): Maximum percentage of read length allowed to not be mapped during filtering (default: 0).
* ``--multimap-threshold=<multimap-threshold>`` (_integer_): Maximum number of multiple mappings allowed during filtering (default: 1).
* ``--reject-multimaps`` (_flag_): If set, any read which multimaps to any species' genome will be rejected and not be assigned to any species.
* ``--best`` (_flag_): Adopt a filtering strategy that provides an excellent balance between sensitivity and specificity. Note that specifying this option overrides the values of the ``--mismatch-threshold``, ``--minmatch-threshold`` and ``--multimap-threshold`` options. In addition, ``--reject-multimaps`` is turned off.
* ``--conservative`` (_flag_): Adopt a filtering strategy where minimising the number of reads mis-assigned to the wrong species takes foremost priority. Note that specifying this option overrides the values of the ``--mismatch-threshold``, ``--minmatch-threshold`` and ``--multimap-threshold options``. In addition, ``--reject-multimaps`` is turned on.
* ``--recall`` (_flag_): Adopt a filtering strategy where sensitivity is prioritised over specificity. Note that specifying this option overrides the values of the ``--mismatch-threshold``, ``--minmatch-threshold`` and ``--multimap-threshold`` options. In addition, ``--reject-multimaps`` is turned off.
* ``--permissive`` (_flag_): Adopt a filtering strategy in which sensitivity is maximised. Note that specifying this option overrides the values of the ``--mismatch-threshold``, ``--minmatch-threshold`` and ``--multimap-threshold`` options. In addition, ``--reject-multimaps`` is turned off.

Performance
-----------

These are optional parameters concerning the running of the pipeline.

* ``-t <num-threads> --num-threads=<num-threads>`` (_integer_): Number of threads to use for parallel processing (default: 1).
* ``--run-separation`` (_flag_): If specified, species separation will be run; otherwise scripts to perform separation will be created but not run. If the option ``--run-separation`` is not specified, a Makefile is written to the given output directory, via which all stages of species separation can be run under the user's control. If ``--run-separation`` is specified, however, the Makefile is both written and executed, and all stages of species separation are performed automatically.
* ``--log-level=<log-level>`` (_text parameter_): Sets the minimum severity level at which log messages will be output (one of "debug", "info", "warning", "error" or "critical").
* ``--delete-intermediate`` (_flag_): If specified, intermediate BAM files (contain raw mapped and sorted reads) will be deleted.
* ``--sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>`` (_text parameter_): Specify the temporary directory to be used by 'sambamba sort' (default: ``/tmp``).

[Next: References](references.md)

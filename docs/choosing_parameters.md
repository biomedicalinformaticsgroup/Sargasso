Choosing parameters
===================

While Sargasso includes a number of "pre-packaged" filtering strategies, each providing a particular balance between sensitivity and specificity, users may be interested to explore choices of filtering parameter in greater depth, allowing them to select parameters that are tailored to their own experimental situation and the sequencing technology used. The ``sargasso_parameter_test`` script is provided to aid this process.

Usage
=====

The ``sargasso_parameter_test`` script runs Sargasso a number of times, on a set of provided samples, using different filtering parameters each time. Filtering statistics are gathered, and informative plots of sensitivity and specificity are produced, allowing users to examine the trade-offs between these two measures. 

n.b. it is intended that ``sargasso_parameter_test`` be run on **single-species** samples, thus allowing incorrectly assigned reads to be easily identified.

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

Script parameters are described below:

* ``<data-type>`` (_text parameter_): One of 'dnaseq' or 'rnaseq' (matches the Sargasso parameter of the same name).
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
* ``<species-info>`` (_text parameter_): Alignment tool for the nth species that the input samples are to be separated into (matches Sargasso parameter).

Output
======

Sargasso filtering will be run for every combination of the values of the parameters ``--mismatch-setting``, ``--minmatch-setting`` and ``--multimap-setting``, though initial mapping will only be run once (or not at all, if the flag ``--skip-init-run`` is set, in which case mapping is assumed to have already been run using this script). The main output of ``sargasso_parameter_test`` are two plots – ``counts.{png|pdf}`` and ``percentage.{png|pdf}``. These show, respectively, the number or percentage of reads which are assigned, marked as ambiguous, or rejected, for each (i) sample, (ii) species into which data is being separated, and (iii) combination of the mismatch, minmatch and multimap parameters.

Examining these graphs can help the user to choose which filtering parameter values will achieve an acceptable assignment of data to the correct species of origin, while minimising misassignment to the wrong species.

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

In this case, we are running the test script using a single sample, whose true species of origin is rat (paths to the raw reads files for the sample are contained in the file ``test_sample.tsv``). Sargasso filtering will be run for every combination of the specified mismatch (0, 2, 4), minmatch (0, 2, 4), and multimap (1) settings – a total of 9 runs – separating the input sample into rat and human reads (where any reads assigned to human will be incorrect assignments). The resulting graphs ``counts.png`` and ``percentage.png`` will be written into the output directory ``~/tmp/sargasso_test``.

[Next: References](references.md)

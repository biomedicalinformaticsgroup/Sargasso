Usage reference
===============

The main ``species_separator`` script has a number of command line options to alter the behaviour and function of read assignment:

    species_separator
        [--log-level=<log-level>]
        [--reads-base-dir=<reads-base-dir>] [--num-threads=<num-threads>]
        [--mismatch-threshold=<mismatch-threshold>]
        [--minmatch-threshold=<minmatch-threshold>]
        [--multimap-threshold=<multimap-threshold>]
        [--reject-multimaps] [--best] [--conservative] [--recall]
        [--run-separation]
        [--delete-intermediate]
        [--star-executable=<star-executable>]
        [--sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>]
        <samples-file> <output-dir>
        (<species> <species-star-info>)
        (<species> <species-star-info>)
        ...

These parameters are listed and explained below for reference, grouped by their functionality.

Core
----

These parameters are required as the base minimum for the execution of the pipeline.

* ``<samples-file>`` (_file path_): TSV file giving paths (relative to ``<reads-base-dir>``) of raw RNA-seq read data files for each sample.
* ``<output-dir>`` (_file path_): Output directory into which the Makefile will be written, and in which species separation will be performed.
* ``<species>`` (_text parameter_): Name of the nth species (note that at least two species must be specified).
* ``<species-star-info>`` (_file paths_): STAR information for the nth species. This parameter should either consist of (i) the path to the STAR index directory for the species _or_ (ii) a comma separated list of two elements, the first of which is the path to a GTF annotation file for the species, and the second is the path to a directory containing genome FASTA files for the species.
    
Mapping
-------

These parameters control the mapping of mixed-species RNA-seq reads to reference genomes.

* ``--reads-base-dir=<reads-base-dir>`` (_file path_): Base directory for raw RNA-seq read data files.

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

Performance
-----------

These are optional parameters concerning the running of the pipeline.

* ``-t <num-threads> --num-threads=<num-threads>`` (_integer_): Number of threads to use for parallel processing (default: 1).
* ``--run-separation`` (_flag_): If specified, species separation will be run; otherwise scripts to perform separation will be created but not run. If the option ``--run-separation`` is not specified, a Makefile is written to the given output directory, via which all stages of species separation can be run under the user's control. If ``--run-separation`` is specified, however, the Makefile is both written and executed, and all stages of species separation are performed automatically.
* ``--log-level=<log-level>`` (_text parameter_): Sets the minimum severity level at which log messages will be output (one of "debug", "info", "warning", "error" or "critical").
* ``--delete-intermediate`` (_flag_): If specified, intermediate BAM files (contain raw mapped and sorted reads) will be deleted.
* ``--star-executable=<star-executable>`` (_text parameter_): Specifies the path to the STAR executable; this can be used to run Sargasso with a particular version of STAR (default ``STAR``).
* ``--sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>`` (_text parameter_): Specify the temporary directory to be used by 'sambamba sort' (default: ``/tmp``).

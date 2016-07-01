Usage reference
===============

The main ``species_separator`` script has a number of command line options to alter the behaviour and functionary of read assignment. These parameters are listed and explained below for reference, grouped by their functionality.

Core
----

These parameters are required as the base minimum for the execution of the pipeline.

``<species-one>``
    *Type: text parameter*

    Name of first species.

``<species-two>``
    *Type: text parameter*

    Name of second species.

``<samples-file>``
    *Type: file path*

    TSV file giving paths (relative to ``<reads-base-dir>``) of raw RNA-seq read data files for each sample.

``<output-dir>``
    *Type: file path*

    Output directory into which the Makefile will be written, and in which species separation will be performed.


Administrative
--------------

These parameters display or create information to help you understand the software package.

``-h --help``
    *Type: flag*

    Outputs the pipeline's help documentation.

``-v --version``
    *Type: flag*

    Outputs the pipeline's version.

``--log-level=<log-level>``
    *Type: text parameter*

    Sets the minimum severity level at which log messages will be output (one of debug, info, warning, error or critical).


Performance
-----------

These are optional parameters concerning the running of the pipeline.

``-t <num-threads> --num-threads=<num-threads>``
    *Type: integer*

    Number of threads to use for parallel processing (default: 1).

``--run-separation``
    *Type: flag*

    If specified, species separation will be run; otherwise scripts to perform separation will be created but not run. If the option "--run-separation" is not specified, a Makefile is written to the given output directory, via which all stages of species separation can be run under the user's control. If "--run-separation" is specified however, the Makefile is both written and executed, and all stages of species separation are performed automatically.


I/O
---

The I/O parameters are used to specify directories for data input and output.

``--reads-base-dir=<reads-base-dir>``
    *Type: file path*

    Base directory for raw RNA-seq read data files.

``--s1-gtf=<species-one-gtf-file>``
    *Type: file path*

    GTF annotation file for first species.

``--s2-gtf=<species-two-gtf-file>``
    *Type: file path*

    GTF annotation file for second species.

``--s1-genome-fasta=<species-one-genome-fasta>``
    *Type: file path*

    Directory containing genome FASTA files for first species.

``--s2-genome-fasta=<species-two-genome-fasta>``
    *Type: file path*

    Directory containing genome FASTA files for second species.

``--s1-index=<species-one-star-index>``
    *Type: file path*

    STAR index directory for first species.

``--s2-index=<species-two-star-index>``
    *Type: file path*

    STAR index directory for second species.


Assignment Criteria & Optimisation
----------------------------------

These parameters are used to specify criteria that affect how reads are assigned to either species, either by choosing individual values for each parameter or by selecting pre-configured assignment profiles.

``--mismatch-threshold=<mismatch-threshold>``
    *Type: float*

    Maximum percentage of read bases allowed to be mismatches against the genome during filtering (default: 0).

``--minmatch-threshold=<minmatch-threshold>``
    *Type: float*

    Maximum percentage of read length allowed to not be mapped during filtering (default: 0).

``--multimap-threshold=<multimap-threshold>``
    *Type: integer*

    Maximum number of multiple mappings allowed during filtering (default: 1).

``--overhang-threshold=<overhang-threshold>``
    *Type: integer*

    If set, specifies the minimum number of bases that are required on either side of an exon boundary for a read mapping to be accepted.

``--reject-multimaps``
    *Type: flag*

    If set, any read which multimaps to either species' genome will be rejected and not be assigned to either species.

``--best``
    *Type: flag*

    Adopt a filtering strategy that provides an excellent balance between sensitivity and specificity. Note that specifying this option overrides the values of the ``--mismatch-threshold``, ``--minmatch-threshold`` and ``--multimap-threshold`` options. In addition, ``--reject-multimaps`` is turned off.

``--conservative``
    *Type: flag*

    Adopt a filtering strategy where minimising the number of reads mis-assigned to the wrong species takes foremost priority. Note that specifying this option overrides the values of the ``--mismatch-threshold``, ``--minmatch-threshold`` and ``--multimap-threshold options``. In addition, ``--reject-multimaps`` is turned on.

``--recall``
    *Type: flag*

    Adopt a filtering strategy where sensitivity is prioritised over specificity. Note that specifying this option overrides the values of the ``--mismatch-threshold``, ``--minmatch-threshold`` and ``--multimap-threshold`` options. In addition, ``--reject-multimaps`` is turned off.

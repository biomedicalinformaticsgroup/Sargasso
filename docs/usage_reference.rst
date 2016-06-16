Usage reference
===============

The pipeline contains a number of criteria and flags that can be specified to alter the behavour and functionality of the read assignment. These parameters will be listed and explained here for reference. They are grouped by their functionality.

Core
----

These parameters are required as the base minimum for the execution of the pipeline, they are the first four params and are entered sequentially in the following order.

<species-one> - Type: Text Parameter
    Name of first species.

<species-two> - Type: Text Parameter
    Name of second species.

<samples-file> - Type: File path
    TSV file giving paths (relative to <reads-base-dir>) of raw RNA-seq read data files for each sample.

<output-dir> - Type: File path
    Output directory into which Makefile will be written, and in which species separation will be performed.


Administrative
--------------

These parameters display or create information to help you understand this software package.

{help_option_spec} - Type: Flag
    Outputs the pipeline's help documentation            

{ver_option_spec} - Type: Flag   
    Outputs the pipeline's version

{log_option_spec} - Type: Flag
    Writes log to file when executing the read assignment


Performance
-----------

The performance parameters are optional parameters concerning the running of the pipeline.

-t <num-threads> --num-threads=<num-threads> - Type: Integer
    Number of threads to use for parallel processing. [default: 1]

--run-separation - Type: Flag
    If specified, species separation will be run; otherwise scripts to perform separation will be created but not run. If the option "--run-separation" is not specified, a Makefile is written to the given output directory, via which all stages of species separation can be run under the user's control. If "--run-separation" is specified however, the Makefile is both written and executed, and all stages of species separation are performed automatically.


I/O
---

The I/O parameters are used to specify directories for data file input and output.

--reads-base-dir=<reads-base-dir> - Type: File path
    Base directory for raw RNA-seq read data files.

--s1-gtf=<species-one-gtf-file> - Type: File path
    GTF annotation file for first species.

--s2-gtf=<species-two-gtf-file> - Type: File path
    GTF annotation file for second species.

--s1-genome-fasta=<species-one-genome-fasta> - Type: File path
    Directory containing genome FASTA files for first species.

--s2-genome-fasta=<species-two-genome-fasta> - Type: File path
    Directory containing genome FASTA files for second species.

--s1-index=<species-one-star-index> - Type: File path
    STAR index directory for first species.

--s2-index=<species-two-star-index> - Type: File path
    STAR index directory for second species.


Assignment Criteria & Optimisation
----------------------------------

These parameters are used to specify criteria that affect how reads are assigned to either species. Either by choosing individual values for each parameter or by selecting pre-configured assignment profiles.

--mismatch-threshold=<mismatch-threshold> - Type: Integer
    Maximum percentage of read bases allowed to be mismatches against the genome during filtering [default: 0].

--minmatch-threshold=<minmatch-threshold> - Type: Integer
    Maximum percentage of read length allowed to not be mapped during filtering [default: 0].

--multimap-threshold=<multimap-threshold> - Type: Integer
    Maximum number of multiple mappings allowed during filtering [default: 1].

--overhang-threshold=<overhang-threshold> - Type: Integer
    If set, allows specification of the minimum number of bases that are allowed on either side of an exon boundary for a read mapping to be accepted

--reject-multimaps - Type: Flag
    If set, any read which multimaps to either species' genome will be rejected and not be assigned to either species.

--best - Type: Flag
    Adopt a filtering strategy that provides an excellent balance between sensitivity and specificity. Note that specifying this option overrides the values of the mismatch-threshold, minmatch-threshold and multimap-threshold options. In addition, reject-multimaps is turned off.

--conservative - Type: Flag
    Adopt a filtering strategy where minimising the number of reads mis-assigned to the wrong species takes foremost priority. Note that specifying this option overrides the values of the mismatch-threshold, minmatch-threshold and multimap-threshold options. In addition, reject-multimaps is turned on.

--recall - Type: Flag
    Adopt a filtering strategy where sensitivity is prioritised over specificity. Note that specifying this option overrides the values of the mismatch-threshold, minmatch-threshold and multimap-threshold options. In addition, reject-multimaps is turned off.



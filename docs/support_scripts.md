Support scripts
===============

In addition to the main ``species_separator`` script, execution of the *Sargasso* pipeline relies on a number of supporting Python and Bash scripts. Their usage patterns are described here; note, however, that in normal usage these scripts need not be executed directly by the user.

build_star_index (Bash)
-----------------------

Usage::

    build_star_index
        <sequence-fasta-files> <gtf-file> <num-threads> <index-dir>

Build ``STAR`` [STAR]_ indices for a species' genome. ``build_star_index`` is called from the species separation Makefile.

Options:

``<sequence-fasta-files>``
    *Type: list of file paths*

    Space-separated list of genome FASTA files.

``<gtf-file>``
    *Type: file path*

    Path to GTF file containing transcript annotations.

``<num-threads>``
    *Type: integer*

    Number of threads to be used for genome generation.

``<index-dir>``
    *Type: file path*

    Path to directory where genome index files will be stored.

collate_raw_reads (Bash)
------------------------

Usage::

    collate_raw_reads
        <samples> <raw-reads-directory> <reads-dir> <reads-type>
        <raw-read-files-1> <raw-read-files-2>

Assemble links to the FASTQ files containing raw RNA-seq reads for each sample. ``collate_raw_reads`` is called from the species separation Makefile.

``<samples>``
    *Type: text parameter*

    Space-separated list of sample names.

``<raw-reads-directory>``
    *Type: file path*

    Base directory for raw RNA-seq read data files.

``<reads-dir>``
    *Type: file path*

    Directory in which links to raw RNA-seq read files will be collated.

``<reads-type>``
    *Type: text parameter*

    Either "single" for single-end reads, or "paired" for paired-end reads.

``<raw-read-files-1>``
    *Type: list of lists of file paths*

    Space-separated list of comma-separated lists of paths to raw RNA-seq read files. Each comma-separated list should correspond to a sample name in the ``<samples>`` parameter, and paths should be given relative to the ``<raw-reads-directory>`` parameter. In the case of paired-end reads, the read files should correspond to the first read of the pair.

``<raw-read-files-2>``
    *Type: list of lists of file paths*

    Space-separated list of comma-separated list of paths to raw RNA-seq read files. Each comma-separated list should correspond to a sample name in the ``<samples>`` parameter, and paths should be given relative to the ``<raw-reads-directory>`` parameter. In the case of paired-end reads, the read files should correspond to the second read of the pair. In the case of single-end reads, this parameter should be omitted.

filter_control (Python)
-----------------------

Usage::

    filter_control
        [--log-level=<log-level>] [--reject-multimaps]
        <block-dir> <output-dir> <sample-name> <species-one> <species-two>
        <mismatch-threshold> <minmatch-threshold> <multimap-threshold> <overhang-threshold>

Takes as input a directory containing pairs of BAM files, each pair being the result of mapping a set of mixed species RNA-seq reads against the two species' genome (in normal operation, all pairs of BAM files will correspond to a single sample, having been split in pieces for efficiency of filtering). Each pair of BAM files is passed to an instance of the script ``filter_sample_reads``, running on a separate thread, which writes filtered read mappings to a pair of species-specific output BAM files. 

``filter_control`` is called by the script ``filter_reads``.

``--log-level=<log-level>``
    *Type: text parameter*

    Sets the minimum severity level at which log messages will be output (one of "debug", "info", "warning", "error" or "criticial").

``--reject-multimaps``
    *Type: flag*

    If set, any read which multimaps to either species' genome will be rejected and not be assigned to either species.

``<block-dir>``
    *Type: file path*

    Directory containing pairs of mapped read BAM files.

``<output-dir>``
    *Type: file path*

    Directory into which species-separated reads will be written.

``<sample-name>``
    *Type: text parameter*

    Name of RNA-seq sample being processed.

``<species-one>``
    *Type: text parameter*

    Name of first species.

``<species-two>``
    *Type: text parameter*

    Name of second species.

``<mismatch-threshold>``
    *Type: float*

    Maximum percentage of read bases allowed to be mismatches against the genome during filtering.

``<minmatch-threshold>``
    *Type: float*

    Maximum percentage of read length allowed to not be mapped during filtering.

``<multimap-threshold>``
    *Type: integer*

    Maximum number of multi-mappings allowed during filtering.

``<overhang-threshold>``
    *Type: integer*

    The minimum number of bases that are allowed on either side of an exon boundary for a read mapping to be accepted.

filter_reads (Bash)
-------------------

Usage::

    filter_reads
        <species-one> <species-two> <samples>
        <input-dir> <output-dir> <num-threads>
        <mismatch-threshold> <minmatch-threshold> <multimap-threshold>
        <overhang-threshold> <reject-multimaps>

For each sample, take the RNA-seq reads mapping to each genome, and assign them to their correct species of origin. ``filter_reads`` is called by the species separation Makefile.

``<species-one>``
    *Type: text parameter*

    Name of species one.

``<species-two>``
    *Type: text parameter*

    Name of species two.

``<samples>``
    *Type: text parameter*

    Space-separated list of sample names.

``<input-dir>``
    *Type: file path*

    Directory containing, for each sample and each species, name-sorted BAM files containing read mappings for that sample's RNA-seq reads to the species' genome reference.

``<output-dir>``
    *Type: file path*

    Directory into which species-separated BAM files are to be written.

``<num-threads>``
    *Type: integer*

    Number of threads to be used during species separation.

``<mismatch-threshold>``
    *Type: float*

    Maximum percentage of read bases allowed to be mismatches against the genome during filtering.

``<minmatch-threshold>``
    *Type: float*

    Maximum percentage of read length allowed to not be mapped during filtering.

``<multimap-threshold>``
    *Type: integer*

    Maximum number of multi-mappings allowed during filtering.

``<overhang-threshold>``
    *Type: integer*

    The minimum number of bases that are allowed on either side of an exon boundary for a read mapping to be accepted.

``<reject-multimaps>``
    *Type: text parameter*

    If set to "--reject-multimaps", any read which multimaps to either species' genome will be rejected and not be assigned to either species.

filter_sample_reads (Python)
----------------------------

Usage::

    filter_sample_reads
        [--log-level=<log-level>] [--reject-multimaps]
        <species-one> <species-one-input-bam> <species-one-output-bam>
        <species-two> <species-two-input-bam> <species-two-output-bam>
        <mismatch-threshold> <minmatch-threshold> <multimap-threshold> <overhang-threshold>

``filter_sample_reads`` takes two BAM files as input, the results of mapping a set of mixed species RNA-seq reads against the two species' genomes, and determines, where possible, from which species each read or read pair originates. Disambiguated read mappings are written to two species-specific output BAM files. Note that the input BAM files *must* be sorted in read order --- failure to ensure input BAM files are correctly sorted will result in erroneous output.

``filter_sample_reads`` is called by the script ``filter_control``.

``--log-level=<log-level>``
    *Type: text parameter*

    Sets the minimum severity level at which log messages will be output (one of "debug", "info", "warning", "error" or "criticial").

``--reject-multimaps``
    *Type: flag*

    If set, any read which multimaps to either species' genome will be rejected and not be assigned to either species.

``<species-one>``
    *Type: text parameter*

    Name of species one.

``<species-one-input-bam>``
    *Type: file path*

    BAM file containing reads mapped against the first species' genome.

``<species-one-output-bam>``
    *Type: file path*

    BAM file to which read mappings assigned to the first species after filtering will be written.

``<species-two>``
    *Type: text parameter*

    Name of species two.

``<species-two-input-bam>``
    *Type: file path*

    BAM file containing reads mapped against the second species' genome.

``<species-two-output-bam>``
    *Type: file path*

    BAM file to which read mappings assigned to the second species after filtering will be written.

``<mismatch-threshold>``
    *Type: float*

    Maximum percentage of read bases allowed to be mismatches against the genome during filtering.

``<minmatch-threshold>``
    *Type: float*

    Maximum percentage of read length allowed to not be mapped during filtering.

``<multimap-threshold>``
    *Type: integer*

    Maximum number of multi-mappings allowed during filtering.

``<overhang-threshold>``
    *Type: integer*

    The minimum number of bases that are allowed on either side of an exon boundary for a read mapping to be accepted.

map_reads (Bash)
----------------

Usage::

    map_reads
        <species> <samples> <star-indices-dir> <num-threads>
        <input-dir> <output-dir> <reads-type>

For each sample, map raw RNA-seq reads to each species' genome. ``map_reads`` is called by the species separation Makefile.

``<species>``
    *Type: text parameter*

    Space-separated list of species names.

``<samples>``
    *Type: text parameter*

    Space-separated list of sample names.

``<star-indices-dir>``
    *Type: file path*

    Directory containing ``STAR`` index directories for each species (or links to index directories).

``<num-threads>``
    *Type: integer*

    Number of threads to be used by ``STAR`` during read mapping.

``<input-dir>``
    *Type: file path*

    Directory containing per-sample directories, each of which contains links to the input raw RNA-seq read files for that sample..

``<output-dir>``
    *Type: file path*

    Directory into which to write BAM files containing read mappings.

``<reads-type>``
    *Type: text parameter*

    Either "single" for single-end reads, or "paired" for paired-end reads.

sort_reads (Bash)
-----------------

Usage::

    sort_reads
        <species> <samples> <num-threads> <input-dir> <output-dir>

For each sample, sort mapped reads for each species into name order. ``sort_reads`` is called by the species separation Makefile.

``<species>``
    *Type: text parameter*

    Space-separated list of species names.

``<samples>``
    *Type: text parameter*

    Space-separated list of sample names.

``<num-threads>``
    *Type: integer*

    Number of threads to be used by ``sambamba`` [Sambamba]_ during read sorting.

``<input-dir>``
    *Type: file path*

    Directory containing BAM files containing read mappings for each sample and species.

``<output-dir>``
    *Type: file path*

    Directory into which to write name-ordered BAM files containing read mappings.

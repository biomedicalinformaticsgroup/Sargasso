Support scripts
===============

In addition to the main ``species_separator`` script, execution of the *Sargasso* pipeline relies on a number of supporting Python and Bash scripts. Their usage patterns are described here; note, however, that in normal usage these scripts need not be executed directly by the user.

build_bowtie2_index (Bash)
--------------------------

Usage:

    build_bowtie2_index
        <sequence-fasta-file> <num-threads> <index-dir> <bowtie2-build-executable>

Build a [Bowtie2](references.md) index for a species' genome. ``build_bowtie2_index`` is called from the species separation Makefile.

Options:

* ``<sequence-fasta-file>`` (_file path_): path to a FASTA file containing genome sequences.
* ``<num-threads>`` (_integer_): Number of threads to be used for genome generation.
* ``<index-dir>`` (_file path_): Path to directory where genome index files will be stored.
* ``<bowtie2-build-executable>`` (_file path_): Path to, or name of, ``bowtie2-build`` executable.

build_star_index (Bash)
-----------------------

Usage:

    build_star_index
        <sequence-fasta-files> <gtf-file> <num-threads> <index-dir> <star-executable>

Build a [STAR](references.md) index for a species' genome. ``build_star_index`` is called from the species separation Makefile.

Options:

* ``<sequence-fasta-files>`` (_list of file paths_): Space-separated list of genome FASTA files.
* ``<gtf-file>`` (_file path_): Path to GTF file containing transcript annotations.
* ``<num-threads>`` (_integer_): Number of threads to be used for genome generation.
* ``<index-dir>`` (_file path_): Path to directory where genome index files will be stored.
* ``<star-executable>`` (_file path_): Path to, or name of, STAR executable.

collate_raw_reads (Bash)
------------------------

Usage:

    collate_raw_reads
        <samples> <raw-reads-directory> <reads-dir> <reads-type>
        <raw-read-files-1> <raw-read-files-2>

Assemble links to the FASTQ files containing raw sequencing reads for each sample. ``collate_raw_reads`` is called from the species separation Makefile.

* ``<samples>`` (_text parameter_): Space-separated list of sample names.
* ``<raw-reads-directory>`` (_file path_): Base directory for raw sequencing read data files.
* ``<reads-dir>`` (_file path_): Directory in which links to raw sequencing read files will be collated.
* ``<reads-type>`` (_text parameter_): Either "single" for single-end reads, or "paired" for paired-end reads.
* ``<raw-read-files-1>`` (_list of lists of file paths_): Space-separated list of comma-separated lists of paths to raw sequencing read files. Each comma-separated list should correspond to a sample name in the ``<samples>`` parameter, and paths should be given relative to the ``<raw-reads-directory>`` parameter. In the case of paired-end reads, the read files should correspond to the first read of the pair.
* ``<raw-read-files-2>`` (_list of lists of file paths_): Space-separated list of comma-separated list of paths to raw sequencing read files. Each comma-separated list should correspond to a sample name in the ``<samples>`` parameter, and paths should be given relative to the ``<raw-reads-directory>`` parameter. In the case of paired-end reads, the read files should correspond to the second read of the pair. In the case of single-end reads, this parameter should be omitted.

filter_control (Python)
-----------------------

Usage:

    filter_control
        [--log-level=<log-level>] [--reject-multimaps]
        <block-dir> <output-dir> <sample-name> 
        <mismatch-threshold> <minmatch-threshold> <multimap-threshold> 
        (<species>) (<species>) ...

Takes as input a directory containing sets of BAM files, each set being the result of mapping a set of mixed species sequencing reads against each species' genome (in normal operation, all pairs of BAM files will correspond to a single sample, having been split in pieces for efficiency of filtering). Each set of BAM files is passed to an instance of the script ``filter_sample_reads``, running on a separate thread, which writes filtered read mappings to a set of species-specific output BAM files. 

``filter_control`` is called by the script ``filter_reads``.

* ``--log-level=<log-level>`` (_text parameter_): Sets the minimum severity level at which log messages will be output (one of "debug", "info", "warning", "error" or "criticial").
* ``--reject-multimaps`` (_flag_): If set, any read which multimaps to either species' genome will be rejected and not be assigned to either species.
* ``<block-dir>`` (_file path_): Directory containing pairs of mapped read BAM files.
* ``<output-dir>`` (_file path_): Directory into which species-separated reads will be written.
* ``<sample-name>`` (_text parameter_): Name of sample being processed.
* ``<mismatch-threshold>`` (_float_): Maximum percentage of read bases allowed to be mismatches against the genome during filtering.
* ``<minmatch-threshold>`` (_float_): Maximum percentage of read length allowed to not be mapped during filtering.
* ``<multimap-threshold>`` (_integer_): Maximum number of multi-mappings allowed during filtering.
* ``<species>`` (_text parameter_): Name of nth species.

filter_reads (Bash)
-------------------

Usage:

    filter_reads
        <data_type> <samples>
        <input-dir> <output-dir> <num-threads>
        <mismatch-threshold> <minmatch-threshold> <multimap-threshold>
        <reject-multimaps>
        (<species>) (<species>) ...

For each sample, take the sequencing reads mapping to each genome, and assign them to their correct species of origin. ``filter_reads`` is called by the species separation Makefile.

* ``<data-type>`` (_text parameter_): One of "dnaseq" or "rnaseq".
* ``<samples>`` (_text parameter_): Space-separated list of sample names.
* ``<input-dir>`` (_file path_): Directory containing, for each sample and each species, name-sorted BAM files containing read mappings for that sample's RNA-seq reads to the species' genome reference.
* ``<output-dir>`` (_file path_): Directory into which species-separated BAM files are to be written.
* ``<num-threads>`` (_integer_): Number of threads to be used during species separation.
* ``<mismatch-threshold>`` (_float_): Maximum percentage of read bases allowed to be mismatches against the genome during filtering.
* ``<minmatch-threshold>`` (_float_): Maximum percentage of read length allowed to not be mapped during filtering.
* ``<multimap-threshold>`` (_integer_): Maximum number of multi-mappings allowed during filtering.
* ``<reject-multimaps>`` (_text parameter_): If set to "--reject-multimaps", any read which multimaps to any species' genome will be rejected and not be assigned to any species.
* ``<log-level>`` (_text parameter_): Sets the minimum severity level at which log messages will be output (one of "debug", "info", "warning", "error" or "criticial").
 
* ``<species>`` (_text parameter_): Name of nth species.

filter_sample_reads (Python)
----------------------------

Usage:

    filter_sample_reads
        [--log-level=<log-level>] [--reject-multimaps]
        <mismatch-threshold> <minmatch-threshold> <multimap-threshold>
        (<species> <species-input-bam> <species-output-bam>)
        (<species> <species-input-bam> <species-output-bam>) ...

``filter_sample_reads`` takes a set of BAM files as input, the results of mapping a set of mixed species sequencing reads against each species' genome, and determines, where possible, from which species each read or read pair originates. Disambiguated read mappings are written to a set of species-specific output BAM files. Note that the input BAM files *must* be sorted in read order (and should contain mappings for the same set of reads) --- failure to ensure input BAM files are correctly sorted will result in erroneous output.

``filter_sample_reads`` is called by the script ``filter_control``.

* ``--log-level=<log-level>`` (_text parameter_): Sets the minimum severity level at which log messages will be output (one of "debug", "info", "warning", "error" or "criticial").
* ``--reject-multimaps`` (_flag_): If set, any read which multimaps to either species' genome will be rejected and not be assigned to either species.
* ``<mismatch-threshold>`` (_float_): Maximum percentage of read bases allowed to be mismatches against the genome during filtering.
* ``<minmatch-threshold>`` (_float_): Maximum percentage of read length allowed to not be mapped during filtering.
* ``<multimap-threshold>`` (_integer_): Maximum number of multi-mappings allowed during filtering.
* ``<species>`` (_text parameter_): Name of nth species.
* ``<species-input-bam>`` (_file path_): BAM file containing reads mapped against the nth species' genome.
* ``<species-output-bam>`` (_file path_): BAM file to which read mappings assigned to the nth species after filtering will be written.

map_reads_dnaseq (Bash)
-----------------------

Usage:

    map_reads_dnaseq
        <species> <samples> <bowtie-indexes-dir> <num-threads>
        <input-dir> <output-dir> <reads-type> <bowtie2-executable>

For each sample, map raw sequencing reads to each species' genome. ``map_reads_dnaseq`` is called by the species separation Makefile.

* ``<species>`` (_text parameter_): Space-separated list of species names.
* ``<samples>`` (_text parameter_): Space-separated list of sample names.
* ``<star-indexes-dir>`` (_file path_): Directory containing Bowtie2 index directories for each species (or links to index directories).
* ``<num-threads>`` (_integer_): Number of threads to be used by Bowtie2 during read mapping.
* ``<input-dir>`` (_file path_): Directory containing per-sample directories, each of which contains links to the input raw sequencing read files for that sample.
* ``<output-dir>`` (_file path_): Directory into which to write BAM files containing read mappings.
* ``<reads-type>`` (_text parameter_): Either "single" for single-end reads, or "paired" for paired-end reads.
* ``<bowtie2-executable>`` (_file path_): Path to, or name of, the Bowtie2 executable.

map_reads_rnaseq (Bash)
-----------------------

Usage:

    map_reads_rnaseq
        <species> <samples> <bowtie-indexes-dir> <num-threads>
        <input-dir> <output-dir> <reads-type> <star-executable>

For each sample, map raw RNA-seq reads to each species' genome. ``map_reads_rnaseq`` is called by the species separation Makefile.

* ``<species>`` (_text parameter_): Space-separated list of species names.
* ``<samples>`` (_text parameter_): Space-separated list of sample names.
* ``<star-indexes-dir>`` (_file path_): Directory containing STAR index directories for each species (or links to index directories).
* ``<num-threads>`` (_integer_): Number of threads to be used by STAR during read mapping.
* ``<input-dir>`` (_file path_): Directory containing per-sample directories, each of which contains links to the input raw sequencing read files for that sample.
* ``<output-dir>`` (_file path_): Directory into which to write BAM files containing read mappings.
* ``<reads-type>`` (_text parameter_): Either "single" for single-end reads, or "paired" for paired-end reads.
* ``<star-executable>`` (_file path_): Path to, or name of, the STAR executable.

sort_reads (Bash)
-----------------

Usage:

    sort_reads
        <species> <samples> <num-threads> <input-dir> <output-dir> <tmp-dir>

For each sample, sort mapped reads for each species into name order. ``sort_reads`` is called by the species separation Makefile.

* ``<species>`` (_text parameter_): Space-separated list of species names.
* ``<samples>`` (_text parameter_): Space-separated list of sample names.
* ``<num-threads>`` (_integer_): Number of threads to be used by ``sambamba`` [Sambamba](references.md) during read sorting.
* ``<input-dir>`` (_file path_): Directory containing BAM files containing read mappings for each sample and species.
* ``<output-dir>`` (_file path_): Directory into which to write name-ordered BAM files containing read mappings.
* ``<tmp-dir>`` (_file path_): Temporary directory to be used by ``sambamba``.

[Next: References](references.md)

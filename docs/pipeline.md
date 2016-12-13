Pipeline description
====================

Separation of mixed-species RNA-seq reads is performed in a number of stages. While in normal usage these steps can be executed with a single command, the *Sargasso* pipeline also affords fine control over execution, should this be desired. Here we describe each stage of species separation in more detail. 

Makefile
--------

The *Sargasso* pipline is invoked through execution of its main Python script, ``species_separator``. This writes a Makefile with targets corresponding to all stages of the pipeline, namely:

* building or linking to ``STAR`` [STAR]_ indices
* collating raw reads files
* mapping reads from all samples to each genome
* sorting mapped reads in preparation for filtering
* filtering mapped reads according to their true species of origin

In this fashion, fine user control over the pipeline is allowed --- processing can be halted and resumed, stages run separately, or particular stages re-run with alterations to parameter values. However, in typical usage, supplying the ``--run-separation`` option to the ``species_separator`` script will cause the the Makefile's main target to be executed immediately after the file has been written.

STAR indices
------------

The *Sargasso* pipeline uses the efficient and accurate ``STAR`` short RNA-seq read aligner to map mixed species reads to reference genomes. For best speed of operation, paths to pre-existing ``STAR`` indices for the two reference genomes can be supplied to the ``species_separator`` script via its ``--s1-index`` and ``--s2-index`` command-line options. Alternatively, if paths to directories containing genome FASTA files for each species are supplied (via the ``--s1-genome`` and ``--s2-genome`` options) and GTF annotation files (via ``--s1-gtf`` and ``-s2-gtf``), then ``STAR`` will be invoked in index-generation mode to build genome indices for each species.

Collating raw reads
-------------------

The path to a TSV file specifying, in turn, the paths to the FASTQ files containing raw RNA-seq reads for each sample being studied should be provided to the ``species_separator`` script through the required ``<samples-file>`` parameter. Checks are made that each raw reads file exists, and links are made to these files within the species separation output directory.

Mapping reads
-------------

Next, the *Sargasso* pipeline maps all reads to each species' genome using the ``STAR`` read aligner. Although we chose to use ``STAR`` for its speed and the accuracy of alignments produced, it would be possible to substitute another mapping tool by altering the read mapping target in the Makefile --- subsequent analysis steps require only BAM files containing the aligned reads.

Note that when invoking ``STAR``, reads are mapped allowing alignments to multiple locations (``--outFilterMultimapNmax 10000``), however only those mappings with an alignment score equal to the maximum are retained (``--outFilterMultimapScoreRange 0``).

Sorting reads
-------------

Mapped RNA-seq reads are subsequently sorted into name order, so that, when filtering according to their true species of origin, the mappings for each read (or each read pair, in the case of paired-end reads) to each genome can be assessed together. Reads are sorted using the ``sambamba`` [Sambamba]_ alignment processing tool.

Filtering reads
---------------

Once reads have been sorted into name order, the mappings to both genomes for each read or read pair are examined in turn, to determine that read or read pair's true species of origin.

Firstly, if a read has alignments to one species' genome, but none to the other, it is provisionally assigned to that species; note, however, that the lack of alignments to a species' genome does not necessarily preclude that species being the read's source. Thus a read's mappings to its provisionally assigned species of origin must satisfy a number of user-defined thresholds:

* *Number of multi-maps*: Controlled by the command-line option ``--multimap-threshold``, the read must have at most this number of alignments to the genome of its putative species of origin (multiple mappings of low quality to one species' genome may indicate that the read's true species of origin is the other, but that its locus in that species' genome reference is missing or incorrectly sequenced).
* *Maximum number of mismatches*: Controlled by the command-line option ``--mismatch-threshold``, none of the read's alignments may exceed this percentage of their length mismatched with respect to the reference genome (a read which maps only poorly to one genome may indicate that the true alignment lies in a missing or poorly defined region of the other species' genome).
* *Minimum length of alignment*: Controlled by the command-line option ``--minmatch-threshold``, all the read's mappings must exceed this percentage of the read's length, including both correct matches and mismatches with respect to the reference genome; that is a minimum percentage length of clipping of the read is tolerated (once more, clipped alignments may mark that the true origin of the read lies in the other species' genome).
* *Structure of alignment*: Controlled by the command-line option ``--overhang-threshold``, for spliced alignments, the minimum length of sequence aligned to any one exon must be at least this number of bases; we have found that alignments with very short exon overhangs were frequently the source of incorrect species assignment.

If all alignments for a read or read pair satisfy these criteria, the read or read pair is assigned to that species; the threshold values can be chosen to balance between the precision of assignment of reads to their correct species of origin, and the recall of the maximum number of reads.

In many cases, however, a read or read pair will align, perhaps in multiple locations, to both species' genomes. In this situation the alignments to each genome are examined and compared. First, each set of alignments is tested against the thresholds listed above. If the criteria are violated for the mappings to both genomes, the read is rejected as it cannot be confidently assigned to either species. If the criteria are violated for the mappings to one species' genome, but satisified for the other, the read is assigned to the latter species.

If the mappings to both species' genomes satisfy all thresholds, these sets of mappings are directly compared:

* If the minimum number of mismatches with respect to the reference genome amongst all mappings in the set is fewer for one species than for the other, the read or read pair is assigned to the former species.
* If the minimum number of mismatches is the same for each set, the length of alignments is examined. If the mappings to one species' genome are of full read length, whereas the other set contains alignments that are clipped, the read or read pair is assigned to the species with full length mappings.
* If all the alignments to both species' genomes are full length, or both sets of mappings contain alignments which are clipped, then the number of multi-mappings of the read or read pair to each species' genome are compared. The genome with fewer multi-mappings is chosen as the likely species of origin.
* Finally, if aligning the read or read pair to each species's genome gives the same number of multi-maps, then the read is discarded as ambiguous.

Different filtering strategies, providing a particular balance between sensitivity and specificity, can be adopted by choice of values of threshold options. While these strategies can be fine-tuned by the user, a number of "pre-packaged" strategies are also available. By specifying the ``--best`` command-line option, a filtering strategy is used that provides an excellent balance between precision and recall in a wide variety of situations, whichver the species of origin of the mixed-species RNA-seq data (note that specifying this, or any of the pre-packaged strategies, overrides the values of the ``--mismatch-threshold``, ``--minmatch-threshold`` and ``--multimap-threshold`` options).

On the other hand, in some cases it may be of particular importance to minimise the number of reads mis-assigned to the wrong species, or to prioritise sensitivity over specificity. These two strategies can be adopted via the ``--conservative`` and ``--recall`` command-line options.

At the end of the filtering stage, a BAM file will have been written for each RNA-seq sample, and for each species, containing the genome alignments of the reads from the sample which were assigned to that species.

Efficiency
----------

In order that the *Sargasso* pipeline operates efficiently, multiple cores can be used wherever possible. Both the ``STAR`` read aligner and the ``sambamba`` alignment processing tool are multi-threaded, and multiple cores are used during species assignment by splitting the input alignment files in chunks and executing filtering in parallel.

The number of cores available at all stages of the pipeline is specified by the ``--num-threads`` command-line option to the ``species_separator`` script.

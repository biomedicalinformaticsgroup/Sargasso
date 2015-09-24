# TODO: At each stage log output to an overall log file

NUM_THREADS=1

SAMPLES=1A1 2A1 3A1 1C1 2C1 3C1 1N1 2N1 3N1
RAW_READ_FILES=/srv/data/ghardingham/neuron_astrocyte/rnaseq/1A1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-001_TP-D5-007_1.sanfastq.gz /srv/data/ghardingham/neuron_astrocyte/rnaseq/2A1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-004_TP-D5-003_1.sanfastq.gz /srv/data/ghardingham/neuron_astrocyte/rnaseq/3A1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-003_TP-D5-002_1.sanfastq.gz /srv/data/ghardingham/neuron_astrocyte/rnaseq/1C1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-001_TP-D5-001_1.sanfastq.gz /srv/data/ghardingham/neuron_astrocyte/rnaseq/2C1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-003_TP-D5-005_1.sanfastq.gz /srv/data/ghardingham/neuron_astrocyte/rnaseq/3C1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-002_TP-D5-002_1.sanfastq.gz /srv/data/ghardingham/neuron_astrocyte/rnaseq/1N1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-001_TP-D5-004_1.sanfastq.gz /srv/data/ghardingham/neuron_astrocyte/rnaseq/2N1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-003_TP-D5-008_1.sanfastq.gz /srv/data/ghardingham/neuron_astrocyte/rnaseq/3N1/150129_M00689_0185_000000000-ABAEY_1_TP-D7-002_TP-D5-005_1.sanfastq.gz

SPECIES_ONE=mouse

# These definitions will only be present if we're building STAR indices from scratch
#SPECIES_ONE_GENOME_FASTA=/srv/data/genome/mouse/ensembl-80/primary_assembly/10.fa
#SPECIES_ONE_GTF=/srv/data/genome/mouse/ensembl-80/Mus_musculus.GRCm38.80.gtf

# This definition will only be present if we're using a pre-built STAR index
SPECIES_ONE_GENOME_DIR=/srv/data/genome/mouse/ensembl-80/STAR_Indices/primary_assembly

SPECIES_TWO=rat

# These definitions will only be present if we're building STAR indices from scratch
#SPECIES_TWO_GENOME_FASTA=/srv/data/genome/rat/ensembl-80/top_level/10.fa
#SPECIES_TWO_GTF=/srv/data/genome/rat/ensembl-80/Rattus_norvegicus.Rnor_6.0.80.gtf

# This definition will only be present if we're using a pre-built STAR index
SPECIES_TWO_GENOME_DIR=/srv/data/genome/rat/ensembl-80/STAR_Indices/top_level

STAR_INDICES=star_indices
STAR_INDEX_MOUSE_MASK=star_index_mouse_mask
STAR_INDEX_RAT_MASK=star_index_rat_mask
COLLATE_RAW_READS=raw_reads
MASKED_READS=masked_reads
MAPPED_READS=mapped_reads
SORTED_READS=sorted_reads
FILTERED_READS=filtered_reads

.PHONY : all clean

all: $(FILTERED_READS)

$(FILTERED_READS): $(SORTED_READS)
	# For each sample, take the reads mapping to each genome and filter them to
	# their correct species of origin
	mkdir -p $(FILTERED_READS)
	filter_reads $(SPECIES_ONE) $(SPECIES_TWO) "$(SAMPLES)" $(SORTED_READS) $(FILTERED_READS)

$(SORTED_READS): $(MAPPED_READS)
	# For each sample, sort the mapped reads into read name order
	mkdir -p $(SORTED_READS)
	sort_reads "$(SPECIES_ONE) $(SPECIES_TWO)" "$(SAMPLES)" $(NUM_THREADS) $(MAPPED_READS) $(SORTED_READS)

$(MAPPED_READS): $(MASKED_READS) $(STAR_INDICES)/$(SPECIES_ONE) $(STAR_INDICES)/$(SPECIES_TWO)
	# TODO: Map 'masked reads' for each sample to each species' genome; mapped
	# reads need to be output sorted by read name
	mkdir -p $(MAPPED_READS)
	map_reads "$(SPECIES_ONE) $(SPECIES_TWO)" "$(SAMPLES)" $(RAW_READS) $(MAPPED_READS)

$(MASKED_READS): $(STAR_INDEX_MOUSE_MASK) $(STAR_INDEX_RAT_MASK) $(COLLATE_RAW_READS)
	# TODO: Map reads in each sample to sequences for each species for which
	# mapping reads should be discarded; retain only those reads which don't
	# map to these sequences in either species. These 'masked' reads should be
	# written to the 'masked_reads' directory.
	mkdir -p $(MASKED_READS)

$(COLLATE_RAW_READS):
	# Create a directory with sub-directories for each sample, each of which
	# contains links to the input raw reads files for that sample
	mkdir -p $(COLLATE_RAW_READS)
	collate_raw_reads "$(SAMPLES)" "$(RAW_READ_FILES)" $(COLLATE_RAW_READS)

$(STAR_INDEX_MOUSE_MASK):
	# TODO: Build STAR index for mouse sequences for which mapping reads should
	# be discarded
	mkdir -p $(STAR_INDEX_MOUSE_MASK)

$(STAR_INDEX_RAT_MASK):
	# TODO: Build STAR index for rat sequences for which mapping reads should
	# be discarded
	mkdir -p $(STAR_INDEX_RAT_MASK)

# This task will only be present if we're building STAR indices from scratch
#$(STAR_INDICES)/$(SPECIES_ONE):
	#mkdir -p $(STAR_INDICES)/$(SPECIES_ONE)
	#build_genome_index $(SPECIES_ONE_GENOME_FASTA) $(SPECIES_ONE_GTF) $(NUM_THREADS) $(STAR_INDICES)/$(SPECIES_ONE)

# This task will only be present if we're using a pre-built STAR index
$(STAR_INDICES)/$(SPECIES_ONE):
	ln -s $(SPECIES_ONE_GENOME_DIR) $(STAR_INDICES)/$(SPECIES_ONE)

# This task will only be present if we're building STAR indices from scratch
#$(STAR_INDICES)/$(SPECIES_TWO):
	#mkdir -p $(STAR_INDICES)/$(SPECIES_TWO)
	#build_genome_index $(SPECIES_TWO_GENOME_FASTA) $(SPECIES_TWO_GTF) $(NUM_THREADS) $(STAR_INDICES)/$(SPECIES_TWO)

# This task will only be present if we're using a pre-built STAR index
$(STAR_INDICES)/$(SPECIES_TWO):
	ln -s $(SPECIES_TWO_GENOME_DIR) $(STAR_INDICES)/$(SPECIES_TWO)

clean:
	rm -rf $(STAR_INDICES)/$(SPECIES_ONE)
	rm -rf $(STAR_INDICES)/$(SPECIES_TWO)
	rm -rf $(STAR_INDEX_MOUSE_MASK)
	rm -rf $(STAR_INDEX_RAT_MASK)
	rm -rf $(COLLATE_RAW_READS)
	rm -rf $(MASKED_READS)
	rm -rf $(MAPPED_READS)
	rm -rf $(SORTED_READS)
	rm -rf $(FILTERED_READS)

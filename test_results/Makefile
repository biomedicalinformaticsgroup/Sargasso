# TODO: At each stage log output to an overall log file

SAMPLES=1A1 2A1 3A1 1C1 2C1 3C1 1N1 2N1 3N1
SPECIES_ONE=mouse
SPECIES_TWO=rat
NUM_THREADS=1

STAR_INDEX_MOUSE_GENOME=star_index_mouse_genome
STAR_INDEX_RAT_GENOME=star_index_rat_genome
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

$(MAPPED_READS): $(MASKED_READS) $(STAR_INDEX_MOUSE_GENOME) $(STAR_INDEX_RAT_GENOME)
	# TODO: Map 'masked reads' for each sample to each species' genome; mapped
	# reads need to be output sorted by read name
	mkdir -p $(MAPPED_READS)
	map_reads "$(SPECIES_ONE) $(SPECIES_TWO)" "$(SAMPLES)" $(MAPPED_READS)

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

$(STAR_INDEX_MOUSE_MASK):
	# TODO: Build STAR index for mouse sequences for which mapping reads should
	# be discarded
	mkdir -p $(STAR_INDEX_MOUSE_MASK)

$(STAR_INDEX_RAT_MASK):
	# TODO: Build STAR index for rat sequences for which mapping reads should
	# be discarded
	mkdir -p $(STAR_INDEX_RAT_MASK)

$(STAR_INDEX_MOUSE_GENOME):
	# TODO: Build STAR index for species 1
	ln -s /srv/data/genome/mouse/ensembl-80/STAR_Indices/primary_assembly $(STAR_INDEX_MOUSE_GENOME)

$(STAR_INDEX_RAT_GENOME):
	# TODO: Build STAR index for species 2
	ln -s /srv/data/genome/rat/ensembl-80/STAR_Indices/top_level $(STAR_INDEX_RAT_GENOME)

clean:
	rm -rf $(STAR_INDEX_MOUSE_GENOME)
	rm -rf $(STAR_INDEX_RAT_GENOME)
	rm -rf $(STAR_INDEX_MOUSE_MASK)
	rm -rf $(STAR_INDEX_RAT_MASK)
	rm -rf $(COLLATE_RAW_READS)
	rm -rf $(MASKED_READS)
	rm -rf $(MAPPED_READS)
	rm -rf $(SORTED_READS)
	rm -rf $(FILTERED_READS)

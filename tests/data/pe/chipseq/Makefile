NUM_THREADS=2

DATA_TYPE=chipseq

MAPPER_EXECUTABLE=bowtie2

SAMBAMBA_SORT_TMP_DIR=/tmp

SAMPLES=chiseq_mouse_sample
RAW_READS_DIRECTORY=/home/xinhe/Projects/Sargasso/tests/data/raw_reads
RAW_READS_FILES_1=chipseq_mouse_pe_R1.fastq.gz
RAW_READS_FILES_2=chipseq_mouse_pe_R2.fastq.gz

MOUSE_BOWTIE2_DIR=/srv/data/genome/mouse/ensembl-94/bowtie2_indices/primary_assembly
HUMAN_BOWTIE2_DIR=/srv/data/genome/human/ensembl-94/bowtie2_indices/primary_assembly
RAT_BOWTIE2_DIR=/srv/data/genome/rat/ensembl-94/bowtie2_indices/toplevel

MAPPER_INDICES=mapper_indexes
COLLATE_RAW_READS=raw_reads
MAPPED_READS=mapped_reads
SORTED_READS=sorted_reads
FILTERED_READS=filtered_reads

.PHONY: all clean

all: $(FILTERED_READS)

$(FILTERED_READS): $(SORTED_READS)
	# For each sample, take the reads mapping to each genome and filter them
	# to their correct species of origin
	mkdir -p $(FILTERED_READS)
	filter_reads $(DATA_TYPE) "$(SAMPLES)" $(SORTED_READS) $(FILTERED_READS) $(NUM_THREADS) 1.0 2.0 999999 "" mouse human rat

$(SORTED_READS): $(MAPPED_READS)
	# For each sample, sort the mapped reads into read name order
	mkdir -p $(SORTED_READS)
	sort_reads "mouse human rat" "$(SAMPLES)" $(NUM_THREADS) $(MAPPED_READS) $(SORTED_READS) $(SAMBAMBA_SORT_TMP_DIR)

$(MAPPED_READS): $(MAPPER_INDICES)/mouse $(MAPPER_INDICES)/human $(MAPPER_INDICES)/rat $(COLLATE_RAW_READS)
	# Map reads for each sample to each species' genome
	mkdir -p $(MAPPED_READS)
	map_reads_chipseq "mouse human rat" "$(SAMPLES)" $(MAPPER_INDICES) $(NUM_THREADS) $(COLLATE_RAW_READS) $(MAPPED_READS) paired bowtie2

$(COLLATE_RAW_READS): 
	# Create a directory with sub-directories for each sample, each of which
	# contains links to the input raw reads files for that sample
	mkdir -p $(COLLATE_RAW_READS)
	collate_raw_reads "$(SAMPLES)" $(RAW_READS_DIRECTORY) $(COLLATE_RAW_READS) paired "$(RAW_READS_FILES_1)" "$(RAW_READS_FILES_2)"

$(MAPPER_INDICES)/mouse: 
	mkdir -p $(MAPPER_INDICES)
	ln -s $(MOUSE_BOWTIE2_DIR) $(MAPPER_INDICES)/mouse

$(MAPPER_INDICES)/human: 
	mkdir -p $(MAPPER_INDICES)
	ln -s $(HUMAN_BOWTIE2_DIR) $(MAPPER_INDICES)/human

$(MAPPER_INDICES)/rat: 
	mkdir -p $(MAPPER_INDICES)
	ln -s $(RAT_BOWTIE2_DIR) $(MAPPER_INDICES)/rat


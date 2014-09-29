REF = ../all_mito2.fas

#This gives all possible accessions, but not all may be downloaded:
ALL_ACCS = $(shell cat ../../raw_data/read_accessions.txt)
#This gives all the possible raw files on disc...
RAW1_ACCS = $(patsubst ../../raw_data/%_1.fastq.gz,%,$(wildcard ../../raw_data/*_1.fastq.gz))
RAW2_ACCS = $(patsubst ../../raw_data/%_2.fastq.gz,%,$(wildcard ../../raw_data/*_2.fastq.gz))
#Take the intersection for accessions with both /1 and /2 reads:
ACCS = $(filter $(filter $(ALL_ACCS), $(RAW1_ACCS)), $(RAW2_ACCS))

#Can now construct the list of bloom-filtered reads, etc
MRFAST_COVERAGE_PNG = $(patsubst %,%_vs_all_mito_mrfast_repair.cov.png,$(ACCS))
MRFAST_MAPPED_READS = $(patsubst %,%_vs_all_mito_mrfast_dedup.bam,$(ACCS))
FILTERED_READS = $(patsubst %,../bloom_filter/%.fastq.bgz,$(ACCS))
INDEXES=$(FILTERED_READS:.fastq.bgz=.fastq.bgz.idx) 

#Define some tools (symlink to picobio github repository?)
INTERLACE_FASTQ = ../../scripts/interlace_fastq.py
BLOOMING_READS = ../../scripts/blooming_reads.py
SEQIO_INDEX = ../../scripts/seqio_index_db.py

all: mrfast_map mrfast_cov
	@echo ACCS: $(ACCS)

#####################################################
mrfast_map: $(MRFAST_MAPPED_READS) $(FILTERED_READS)

%_vs_all_mito_mrfast_dedup.bam: ../bloom_filter/%.fastq.bgz $(REF)
	@echo "Mapping filtering reads for $* to make $@"

#####################################################

mrfast_cov: $(MRFAST_COVERAGE_PNG) $(MAPPED_READS)

%_vs_all_mito_mrfast_repair.cov.png: %_vs_all_mito_mrfast_dedup.bam
	@echo "Plotting coverage for $* to make $@"

REF = ../all_mito2.fas

#This gives all possible accessions, but not all may be downloaded:
ALL_ACCS = $(shell cat ../../raw_data/read_accessions.txt)
#This gives all the possible raw files on disc...
RAW1_ACCS = $(patsubst ../../raw_data/%_1.fastq.gz,%,$(wildcard ../../raw_data/*_1.fastq.gz))
RAW2_ACCS = $(patsubst ../../raw_data/%_2.fastq.gz,%,$(wildcard ../../raw_data/*_2.fastq.gz))
#Take the intersection for accessions with both /1 and /2 reads:
ACCS = $(filter $(filter $(ALL_ACCS), $(RAW1_ACCS)), $(RAW2_ACCS))

#Can now construct the list of bloom-filtered reads, etc
RAW1 = $(patsubst %,../../raw_data/%_1.fastq.gz,$(ACCS))
RAW2 = $(patsubst %,../../raw_data/%_2.fastq.gz,$(ACCS))
FILTERED_READS = $(patsubst %,%.fastq.bgz,$(ACCS))
INDEXES=$(FILTERED_READS:.fastq.bgz=.fastq.bgz.idx) 

#Define some tools (symlink to picobio github repository?)
INTERLACE_FASTQ = ../../scripts/interlace_fastq.py
BLOOMING_READS = ../../scripts/blooming_reads.py
SEQIO_INDEX = ../../scripts/seqio_index_db.py

all: bloom_filter index_bgz
	@echo ACCS: $(ACCS)

#####################################################
bloom_filter: $(FILTERED_READS)

%.fastq.bgz: ../../raw_data/%_1.fastq.gz ../../raw_data/%_2.fastq.gz
	@echo Filtering raw reads for $* to make $@
	$(INTERLACE_FASTQ) $(patsubst %,../../raw_data/%_1.fastq.gz,$*) $(patsubst %,../../raw_data/%_2.fastq.gz,$*) | $(BLOOMING_READS) -c $(REF) -k 20 -m 1 -f fastq | bgzip > $@

#####################################################

index_bgz: $(INDEXES)

%.fastq.bgz.idx: %.fastq.bgz
	@echo (Re)indexing filtered reads for  $* to make $@
	rm -f $@
	$(SEQIO_INDEX) -f fastq $^

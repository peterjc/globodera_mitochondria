REF = ../all_mito2.fas

#This gives all possible accessions, but not all may be downloaded:
ALL_ACCS = $(shell cat ../../raw_data/read_accessions.txt)
#This gives all the possible raw files on disc...
RAW1_ACCS = $(patsubst ../../raw_data/%_1.fastq.gz,%,$(wildcard ../../raw_data/*_1.fastq.gz))
RAW2_ACCS = $(patsubst ../../raw_data/%_2.fastq.gz,%,$(wildcard ../../raw_data/*_2.fastq.gz))
#Take the intersection for accessions with both /1 and /2 reads:
ACCS = $(filter $(filter $(ALL_ACCS), $(RAW1_ACCS)), $(RAW2_ACCS))

#Can now construct the lists for paired reads
P1 = $(patsubst %,../trimmed/%_1P.fastq.gz,$(ACCS))
P2 = $(patsubst %,../trimmed/%_2P.fastq.gz,$(ACCS))
FILTERED_P_READS = $(patsubst %,%_P.fastq.bgz,$(ACCS))
P_INDEXES=$(FILTERED_P_READS:.fastq.bgz=.fastq.bgz.idx)

#Lists for orphan reads
U1 = $(patsubst %,../trimmed/%_1U.fastq.gz,$(ACCS))
U2 = $(patsubst %,../trimmed/%_2U.fastq.gz,$(ACCS))
FILTERED_U_READS = $(patsubst %,%_U.fastq.bgz,$(ACCS))

#Define some tools (symlink to picobio github repository?)
INTERLACE_FASTQ = ~/repositories/picobio/blooming_reads/interlace_fastq.py
BLOOMING_READS = ~/repositories/picobio/blooming_reads/blooming_reads.py
SEQIO_INDEX = ~/repositories/picobio/seq_manipulation/seqio_index_db.py

all: bloom_filter index_bgz
	@echo ACCS: $(ACCS)

#####################################################
bloom_filter: $(FILTERED_P_READS) $(FILTERED_U_READS)

$(FILTERED_P_READS): %_P.fastq.bgz: ../trimmed/%_1P.fastq.gz ../trimmed/%_2P.fastq.gz
	@echo "Filtering raw paired reads for $* to make $@"
	$(INTERLACE_FASTQ) $(patsubst %,../trimmed/%_1P.fastq.gz,$*) $(patsubst %,../trimmed/%_2P.fastq.gz,$*) | $(BLOOMING_READS) -p -c $(REF) -k 20 -m 1 -f fastq | bgzip > $@

$(FILTERED_U_READS): %_U.fastq.bgz: ../trimmed/%_1U.fastq.gz ../trimmed/%_2U.fastq.gz
	@echo "Filtering raw orphaned reads for $* to make $@"
	zcat ../trimmed/$*_1U.fastq.gz ../trimmed/$*_2U.fastq.gz | $(BLOOMING_READS) -s -c $(REF) -k 20 -m 1 -f fastq | bgzip > $@

#####################################################

index_bgz: $(P_INDEXES)

$(P_INDEXES): %.fastq.bgz.idx: %.fastq.bgz
	@echo "(Re)indexing filtered reads for $* to make $@"
	rm -f $@
	$(SEQIO_INDEX) -f fastq $^

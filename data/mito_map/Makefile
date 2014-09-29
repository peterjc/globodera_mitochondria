REF1 = ../all_mito.fas
REF1_INDEX = ../all_mito.fas.fai
REF2 = ../all_mito2.fas
REF2_INDEX = ../all_mito2.fas.fai

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
INDEXES = $(FILTERED_READS:.fastq.bgz=.fastq.bgz.idx) 

#Define some tools (symlink to picobio github repository?)
MRFAST = ~/bin/mrfast-2.5.0.0
DEDUP_CIRCLE = ~/repositories/picobio/blooming_reads/dedup_circular_sam.py
REPAIR_CIRCLE = ~/repositories/picobio/blooming_reads/re_pair_circular_sam.py

all: mrfast_map mrfast_cov
	@echo ACCS: $(ACCS)

$(REF1_INDEX): $(REF1)
	@echo "Indexing $(REF1)"
	samtools faidx $@

$(REF2_INDEX): $(REF2)
	@echo "Indexing $(REF2)"
	samtools faidx $@

#####################################################
mrfast_map: $(MRFAST_MAPPED_READS) $(FILTERED_READS) $(REF2_INDEX)

#Intermediate file, raw SAM output for mapping against doubled reference
%_vs_all_mito2_mrfast.sam: ../bloom_filter/%.fastq.bgz $(REF2) $(REF2_INDEX)
	@echo "Mapping $* to doubled-reference using mrfast to make $@"
	$(MRFAST) --seq ../bloom_filter/$*.fastq.bgz --seqcomp --search $(REF2) -o $@

#Intermediate file, raw SAM output deduplicated
%_vs_all_mito_mrfast_dedup_unsorted.bam: %_vs_all_mito2_mrfast.sam $(REF1)
	@echo "De-duplicating $* mapping to make $@"
	$(DEDUP_CIRCLE) -i $*_vs_all_mito2_mrfast.sam -c $(REF1) | samtools view -u -S - > $@

#Output file, raw SAM output deduplicated and sorted
%_vs_all_mito_mrfast_dedup.bam: %_vs_all_mito_mrfast_dedup_unsorted.bam
	@echo "Sorting de-duplicated $* mapping to make $@"
	samtools sort %_vs_all_mito_mrfast_dedup_unsorted.bam $*_vs_all_mito_mrfast_dedup

#####################################################

mrfast_cov: $(MRFAST_COVERAGE_PNG) $(MAPPED_READS)

%_vs_all_mito_mrfast_repair.cov.png: %_vs_all_mito_mrfast_dedup.bam
	@echo "Plotting coverage for $* to make $@"

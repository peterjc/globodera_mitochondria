REF1 = ../all_mito.fas
REF1_INDEX = ../all_mito.fas.bwt
REF2 = ../all_mito2.fas
REF2_INDEX = ../all_mito2.fas.bwt

#This gives all possible accessions, but not all may be downloaded:
ALL_ACCS = $(shell cat ../../raw_data/read_accessions.txt)
#This gives all the possible raw files on disc...
RAW1_ACCS = $(patsubst ../../raw_data/%_1.fastq.gz,%,$(wildcard ../../raw_data/*_1.fastq.gz))
RAW2_ACCS = $(patsubst ../../raw_data/%_2.fastq.gz,%,$(wildcard ../../raw_data/*_2.fastq.gz))
#Take the intersection for accessions with both /1 and /2 reads:
ACCS = $(filter $(filter $(ALL_ACCS), $(RAW1_ACCS)), $(RAW2_ACCS))

#Can now construct the lists for paired reads
P = $(patsubst %,../bloom_filter/%_P.fastq.bgz,$(ACCS))
#Lists for orphan reads
U = $(patsubst %,../bloom_filter/%_U.fastq.bgz,$(ACCS))

BWA_P = $(patsubst %,%_P.sam,$(ACCS))
BWA_U =	$(patsubst %,%_U.sam,$(ACCS))

#Define some tools (symlink to picobio github repository?)
BWA = ~/Downloads/bwa-0.7.10/bwa
DEDUP_CIRCLE = ~/repositories/picobio/blooming_reads/dedup_circular_sam.py
REPAIR_CIRCLE = ~/repositories/picobio/blooming_reads/re_pair_circular_sam.py

all: bwa_map
	@echo ACCS: $(ACCS)

$(REF1_INDEX): $(REF1)
	@echo "BWA indexing $(REF1)"
	bwa index $@

$(REF2_INDEX): $(REF2)
	@echo "BWA indexing $(REF2)"
	bwa index $@

bwa_map: bwa_p_map bwa_u_map

#####################################################
bwa_p_map: $(BWA_P)

bwa_u_map: $(BWA_U)

#Intermediate file, raw SAM output for mapping against doubled reference
$(BWA_P): %_P.sam: ../bloom_filter/%_P.fastq.bgz $(REF2) $(REF2_INDEX)
	@echo "BWA mapping $* paired reads to make $@"
	$(BWA) mem -p -S -a $(REF2) ../bloom_filter/$*_P.fastq.bgz > $*_P.sam

$(BWA_U): %_U.sam: ../bloom_filter/%_U.fastq.bgz $(REF2) $(REF2_INDEX)
	@echo "BWA mapping $* orphaned reads to make $@"
	$(BWA) mem -S -a $(REF2) ../bloom_filter/$*_U.fastq.bgz > $*_U.sam

#####################################################
#
# Example, trying BWA-MEM since it supports interleaved FASTQ
#
# $ time bwa mem -p -S -a ../all_mito2.fas ../bloom_filter/ERR019316_P.fastq.bgz > ERR019316_P.sam
# [M::main_mem] read 148262 sequences (10000090 bp)...
# ...
# [M::main_mem] read 36806 sequences (2542231 bp)...
# [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (0, 0, 0, 0)
# [M::mem_pestat] skip orientation FF as there are not enough pairs
# [M::mem_pestat] skip orientation FR as there are not enough pairs
# [M::mem_pestat] skip orientation RF as there are not enough pairs
# [M::mem_pestat] skip orientation RR as there are not enough pairs
# [M::mem_process_seqs] Processed 36806 reads in 2.988 CPU sec, 2.994 real sec
# [main] Version: 0.7.10-r789
# [main] CMD: bwa mem -p -S -a ../all_mito2.fas ../bloom_filter/ERR019316_P.fastq.bgz
# [main] Real time: 76.536 sec; CPU: 73.088 sec
# real	1m16.555s
# user	1m11.987s
# sys	0m1.117s
#
# And for the unpaired reads,
#
# $ time bwa mem -S -a ../all_mito2.fas ../bloom_filter/ERR019316_U.fastq.bgz  > ERR019316_U.sam
# [M::main_mem] read 145816 sequences (10000139 bp)...
# [M::mem_process_seqs] Processed 145816 reads in 12.073 CPU sec, 12.071 real sec
# [M::main_mem] read 18581 sequences (1235911 bp)...
# [M::mem_process_seqs] Processed 18581 reads in 2.236 CPU sec, 2.235 real sec
# [main] Version: 0.7.10-r789
# [main] CMD: bwa mem -S -a ../all_mito2.fas ../bloom_filter/ERR019316_U.fastq.bgz
# [main] Real time: 15.654 sec; CPU: 14.927 sec
# real	0m15.697s
# user	0m14.700s
# sys	0m0.243s
#
# Example, trying BWA-backtrack (aln/samse/sampe) since reads are up to 100bp
#
# $ bwa aln ../all_mito2.fas ../bloom_filter/ERR019316_P.fastq.bgz -f ERR019316_P.sai
# [bwa_aln] 17bp reads: max_diff = 2
# [bwa_aln] 38bp reads: max_diff = 3
# [bwa_aln] 64bp reads: max_diff = 4
# [bwa_aln] 93bp reads: max_diff = 5
# [bwa_aln] 124bp reads: max_diff = 6
# [bwa_aln] 157bp reads: max_diff = 7
# [bwa_aln] 190bp reads: max_diff = 8
# [bwa_aln] 225bp reads: max_diff = 9
# [bwa_aln_core] calculate SA coordinate... 6.19 sec
# [bwa_aln_core] write to the disk... 0.02 sec
# [bwa_aln_core] 262144 sequences have been processed.
# [bwa_aln_core] calculate SA coordinate... 6.19 sec
# [bwa_aln_core] write to the disk... 0.01 sec
# [bwa_aln_core] 524288 sequences have been processed.
# [bwa_aln_core] calculate SA coordinate... 6.27 sec
# [bwa_aln_core] write to the disk... 0.01 sec
# [bwa_aln_core] 786432 sequences have been processed.
# [bwa_aln_core] calculate SA coordinate... 3.10 sec
# [bwa_aln_core] write to the disk... 0.01 sec
# [bwa_aln_core] 920054 sequences have been processed.
# [main] Version: 0.7.10-r789
# [main] CMD: bwa aln -f ERR019316_P.sai ../all_mito2.fas ../bloom_filter/ERR019316_P.fastq.bgz
# [main] Real time: 24.891 sec; CPU: 24.789 sec
#
# $ bwa samse ../all_mito2.fas ERR019316_P.sai ../bloom_filter/ERR019316_P.fastq.bgz -f ERR019316_P.sam
# [bwa_aln_core] convert to sequence coordinate... 0.65 sec
# [bwa_aln_core] refine gapped alignments... 0.47 sec
# [bwa_aln_core] print alignments... 0.73 sec
# [bwa_aln_core] 262144 sequences have been processed.
# [bwa_aln_core] convert to sequence coordinate... 0.64 sec
# [bwa_aln_core] refine gapped alignments... 0.51 sec
# [bwa_aln_core] print alignments... 0.58 sec
# [bwa_aln_core] 524288 sequences have been processed.
# [bwa_aln_core] convert to sequence coordinate... 0.69 sec
# [bwa_aln_core] refine gapped alignments... 0.46 sec
# [bwa_aln_core] print alignments... 0.57 sec
# [bwa_aln_core] 786432 sequences have been processed.
# [bwa_aln_core] convert to sequence coordinate... 0.32 sec
# [bwa_aln_core] refine gapped alignments... 0.22 sec
# [bwa_aln_core] print alignments... 0.30 sec
# [bwa_aln_core] 920054 sequences have been processed.
# [main] Version: 0.7.10-r789
# [main] CMD: bwa samse -f ERR019316_P.sam ../all_mito2.fas ERR019316_P.sai ../bloom_filter/ERR019316_P.fastq.bgz
# [main] Real time: 10.222 sec; CPU: 9.494 sec
#
# Treated as single end reads, BUT, our /1 and /2 suffices were lost :(

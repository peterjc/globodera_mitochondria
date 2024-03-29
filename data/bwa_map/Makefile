#TODO - Run BAM EOF check

REF1 = ../all_mito.fas
REF1_INDEX = ../all_mito.fas.bwt
REF2 = ../all_mito2.fas
REF2_INDEX = ../all_mito2.fas.bwt

REFS = $(shell cat ../../raw_data/references.txt)
#Uncomment for testing:
#REFS = AJ249395

SINGLE_REF = $(patsubst %,%.fas,$(REFS))
DOUBLE_REF = $(patsubst %,%_X2.fas,$(REFS))

SINGLE_BWT = $(patsubst %,%.fas.bwt,$(REFS))
DOUBLE_BWT = $(patsubst %,%_X2.fas.bwt,$(REFS))

#This gives all possible accessions, but not all may be downloaded:
ALL_ACCS = $(shell cat ../../raw_data/read_accessions.txt)
#Uncomment for testing:
#ALL_ACCS = ERR019315
#This gives all the possible raw files on disc...
RAW1_ACCS = $(patsubst ../../raw_data/%_1.fastq.gz,%,$(wildcard ../../raw_data/*_1.fastq.gz))
RAW2_ACCS = $(patsubst ../../raw_data/%_2.fastq.gz,%,$(wildcard ../../raw_data/*_2.fastq.gz))
#Take the intersection for accessions with both /1 and /2 reads:
ACCS = $(filter $(filter $(ALL_ACCS), $(RAW1_ACCS)), $(RAW2_ACCS))

#Can now construct the lists for paired reads
P = $(patsubst %,../bloom_filter/%_P.fastq.bgz,$(ACCS))
#Lists for orphan reads
U = $(patsubst %,../bloom_filter/%_U.fastq.bgz,$(ACCS))

COMBO = $(foreach a,$(ACCS),$(foreach r,$(REFS),$(a)_vs_$(r)))

get_acc = $(firstword $(subst _vs_, , $1))
get_ref = $(word 2,$(subst _vs_, , $1))

#Want these to be auto-deleted intermediate files...
BWA_X2_P_SAM = $(patsubst %,%_X2_P.sam,$(COMBO))
BWA_X2_U_SAM = $(patsubst %,%_X2_U.sam,$(COMBO))
BWA_X2_P_BAM = $(patsubst %,%_X2_P.bam,$(COMBO))
BWA_X2_U_BAM = $(patsubst %,%_X2_U.bam,$(COMBO))

BWA_X2_BAM = $(patsubst %,%_X2.bam,$(COMBO))
BWA_X2_COV = $(patsubst %,%_X2.cov,$(COMBO))
BWA_X2_COV_PNG = $(patsubst %,%_X2_bwa_coverage.png,$(ACCS))

BWA_BAM = $(patsubst %,%.bam,$(COMBO))
BWA_COV = $(patsubst %,%.cov,$(COMBO))
BWA_COV_PNG = $(patsubst %,%_bwa_coverage.png,$(ACCS))

#Define some tools (symlink to picobio github repository?)
BWA = ~/Downloads/bwa-0.7.10/bwa
CHECK_EOF = ~/repositories/picobio/sambam/bgzf_check_eof.py
DEDUP_CIRCLE = ~/repositories/picobio/sambam/sam_undouble_circles.py
CALC_COVERAGE = ~/repositories/picobio/sambam/sam_circular_coverage.py
COVERAGE_PLOT = ../../scripts/stack_coverage_plot.py

empty :=
space := $(empty) $(empty)

all: single double bwa_map bwa_coverage bwa_plots
	@echo ACCS: $(ACCS)
	@echo REFS: $(REFS)

#These are intermediate, but I will remove them once not needed
.SECONDARY: $(BWA_X2_P_SAM) $(BWA_X2_U_SAM) $(BWA_X2_P_BAM) $(BWA_X2_U_BAM)

single: $(SINGLE_REF) $(SINGLE_BWT)

double: $(DOUBLE_REF) $(DOUBLE_BWT)

$(SINGLE_REF): %.fas: ../../raw_data/%.fasta
	@echo "Missing $@ single reference"
	ln -s ../../raw_data/$*.fasta $*.fas

$(DOUBLE_REF): %_X2.fas: %.fas
	@echo "Missing $@ double reference"
	python -c "from Bio import SeqIO; SeqIO.write((r+r for r in SeqIO.parse('$*.fas', 'fasta')), '$*_X2.fas', 'fasta')"

$(SINGLE_BWT): %.fas.bwt: %.fas
	@echo "BWA indexing $* single reference"
	bwa index $*.fas

$(DOUBLE_BWT): %_X2.fas.bwt: %_X2.fas
	@echo "BWA indexing $* double reference"
	bwa index $*_X2.fas

$(REF1_INDEX): $(REF1)
	@echo "BWA indexing $(REF1)"
	bwa index $@

$(REF2_INDEX): $(REF2)
	@echo "BWA indexing $(REF2)"
	bwa index $@

#####################################################
#For now explicitly keep the X2 BAM files...
bwa_map: $(BWA_BAM)

#intermediate files, raw SAM output from mapping against doubled reference
#TODO - how to set dependency to only the specific files needed?
$(BWA_X2_P_SAM): %_X2_P.sam: $(DOUBLE_REF) $(P) $(U)
	@echo "Mapping $(call get_acc,$*) pairs onto $(call get_ref,$*)"
	$(BWA) mem -p -S -a ../refs/$(call get_ref,$*)_X2.fas ../bloom_filter/$(call get_acc,$*)_P.fastq.bgz > $@ 2> $*_X2_P.sam.log
	rm $*_X2_P.sam.log

#intermediate files, raw SAM output from mapping against doubled reference
#TODO - how to set dependency to only the specific files needed?
$(BWA_X2_U_SAM): %_X2_U.sam: $(DOUBLE_REF) $(P) $(U)
	@echo "Mapping $(call get_acc,$*) orphans onto $(call get_ref,$*)"
	$(BWA) mem -S -a ../refs/$(call get_ref,$*)_X2.fas ../bloom_filter/$(call get_acc,$*)_U.fastq.bgz > $@ 2> $*_X2_U.sam.log
	rm $*_X2_U.sam.log

#intermediate files, raw BAM output from mapping against doubled reference
#Should be able to let make delete the input file at the end of the run, but can do now.
$(BWA_X2_P_BAM): %_X2_P.bam: %_X2_P.sam
	@echo "Converting $< into $@"
	samtools view -S -b $< > $@
	$(CHECK_EOF) $@
	rm $<

#intermediate files, raw BAM output from mapping against doubled reference
#Should be able to let make delete the input file at the end of the run, but can do now.
$(BWA_X2_U_BAM): %_X2_U.bam: %_X2_U.sam
	@echo "Converting $< into $@"
	samtools view -S -b $< > $@
	$(CHECK_EOF) $@
	rm $<

#Semi-intermediate files, merged BAM output from mapping pairs & singletons to doubled ref
#Should be able to let make delete its input files at end of run, but can do now.
$(BWA_X2_BAM): %_X2.bam: %_X2_P.bam %_X2_U.bam
	@echo "Merging $^ into $@"
	samtools merge -f $@ $^
	$(CHECK_EOF) $@
	rm $^

#TODO - Sort at this point?
$(BWA_BAM): %.bam: %_X2.bam
	@echo "Applying de-circle logic and dropping unmapped reads to $< to make $@"
	samtools view -F 4 -h $< | $(DEDUP_CIRCLE) -c ../refs/$(call get_ref,$*).fas | samtools view -S -b - > $@
	$(CHECK_EOF) $@

#####################################################

bwa_coverage: $(BWA_COV) $(BWA_X2_COV)

$(BWA_COV): %.cov: %.bam
	@echo "Computing coverage of $(call get_acc,$*) mapped onto $(call get_ref,$*)"
	$(CHECK_EOF) $<
	samtools view -h $< | $(CALC_COVERAGE) -c ../refs/$(call get_ref,$*).fas -o $@
	@touch $@

$(BWA_X2_COV): %_X2.cov: %_X2.bam
	@echo "Computing coverage of $(call get_acc,$*) mapped onto $(call get_ref,$*) (doubled)"
	$(CHECK_EOF) $<
	samtools view -h $< | $(CALC_COVERAGE) -c ../refs/$(call get_ref,$*).fas -o $@
	@touch $@

bwa_plots: $(BWA_COV_PNG) $(BWA_X2_COV_PNG)

#This claims to depend on *all* to coverage files, actually only
#depends on those for the current read accession. How to say that?
$(BWA_COV_PNG): %_bwa_coverage.png: $(BWA_COV)
	@echo "Plotting $* coverage to make $@"
	$(COVERAGE_PLOT) -o $@ $(foreach r,$(REFS),$*_vs_$(r).cov)
	@touch $@

$(BWA_X2_COV_PNG): %_X2_bwa_coverage.png: $(BWA_X2_COV)
	@echo "Plotting $* coverage to make $@"
	$(COVERAGE_PLOT) -o $@ $(foreach r,$(REFS),$*_vs_$(r)_X2.cov)
	@touch $@

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

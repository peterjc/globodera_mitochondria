#This gives all possible accessions, but not all may be downloaded:
ALL_ACCS = $(shell cat ../../raw_data/read_accessions.txt)
#This gives all the possible raw files on disc...
RAW1_ACCS = $(patsubst ../../raw_data/%_1.fastq.gz,%,$(wildcard ../../raw_data/*_1.fastq.gz))
RAW2_ACCS = $(patsubst ../../raw_data/%_2.fastq.gz,%,$(wildcard ../../raw_data/*_2.fastq.gz))
#Take the intersection for accessions with both /1 and /2 reads:
ACCS = $(filter $(filter $(ALL_ACCS), $(RAW1_ACCS)), $(RAW2_ACCS))

#Can now construct the list of reads, etc
RAW1 = $(patsubst %,../../raw_data/%_1.fastq.gz,$(ACCS))
RAW2 = $(patsubst %,../../raw_data/%_2.fastq.gz,$(ACCS))

TRIM_P1 =  $(patsubst %,%_1P.fastq.gz,$(ACCS))
TRIM_U1 =  $(patsubst %,%_1U.fastq.gz,$(ACCS))
TRIM_P2 =  $(patsubst %,%_2P.fastq.gz,$(ACCS))
TRIM_U2 =  $(patsubst %,%_2U.fastq.gz,$(ACCS))

#Tools etc
TRIMMOMATIC = /opt/Trimmomatic-0.32/trimmomatic-0.32.jar
#Should we look for Illumina adapters? Which ones? Does it vary over read sets?
#... ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10
TRIM_ARGS = LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


all: trim
	@echo ACCS: $(ACCS)

#####################################################
trim: $(TRIM_P1) $(TRIM_U1) $(TRIM_P2) $(TRIM_U2)
bloom_filter: $(FILTERED_READS)

%_1P.fastq.gz: ../../raw_data/%_1.fastq.gz ../../raw_data/%_2.fastq.gz
	@echo "Trimming raw reads for $*"
	java -jar $(TRIMMOMATIC) PE -threads 8 -phred33 ../../raw_data/$*_1.fastq.gz ../../raw_data/$*_2.fastq.gz $*_1P.fastq.gz $*_1U.fastq.gz $*_2P.fastq.gz $*_2U.fastq.gz $(TRIM_ARGS)

trim: $(TRIM_P1) $(TRIM_U1) $(TRIM_P2) $(TRIM_U2)
	@echo "Trim reads"

# Example of a sane looking library:
#
# $ time java -jar /opt/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 ../../raw_data/ERR019317_1.fastq.gz ../../raw_data/ERR019317_2.fastq.gz ERR019317_1P.fastq.gz ERR019317_1U.fastq.gz ERR019317_2P.fastq.gz ERR019317_2U.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# TrimmomaticPE: Started with arguments: -threads 8 -phred33 ../../raw_data/ERR019317_1.fastq.gz ../../raw_data/ERR019317_2.fastq.gz ERR019317_1P.fastq.gz ERR019317_1U.fastq.gz ERR019317_2P.fastq.gz ERR019317_2U.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# Input Read Pairs: 12794791 Both Surviving: 7417590 (57.97%) Forward Only Surviving: 2086891 (16.31%) Reverse Only Surviving: 1852859 (14.48%) Dropped: 1437451 (11.23%)
# TrimmomaticPE: Completed successfully
#
# real	8m24.283s
# user	19m27.332s
# sys	1m15.081s
#
# $ ls -l -h ERR019317*
# -rw-r--r-- 1 pc40583 ics 474M Oct  2 10:32 ERR019317_1P.fastq.gz
# -rw-r--r-- 1 pc40583 ics 129M Oct  2 10:32 ERR019317_1U.fastq.gz
# -rw-r--r-- 1 pc40583 ics 485M Oct  2 10:32 ERR019317_2P.fastq.gz
# -rw-r--r-- 1 pc40583 ics 116M Oct  2 10:32 ERR019317_2U.fastq.gz
#
# Example (which seems to a be poor quality library, most reads left in 1U, lots of CCCCCCCCC reads too!)
#
# $ time java -jar /opt/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 ../../raw_data/ERR019315_1.fastq.gz ../../raw_data/ERR019315_2.fastq.gz ERR019315_1P.fastq.gz ERR019315_1U.fastq.gz ERR019315_2P.fastq.gz ERR019315_2U.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# TrimmomaticPE: Started with arguments: -threads 8 -phred33 ../../raw_data/ERR019315_1.fastq.gz ../../raw_data/ERR019315_2.fastq.gz ERR019315_1P.fastq.gz ERR019315_1U.fastq.gz ERR019315_2P.fastq.gz ERR019315_2U.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# Input Read Pairs: 15672763 Both Surviving: 118 (0.00%) Forward Only Surviving: 13821597 (88.19%) Reverse Only Surviving: 20 (0.00%) Dropped: 1851028 (11.81%)
# TrimmomaticPE: Completed successfully
#
# real	13m7.002s
# user	14m16.748s
# sys	0m43.299s
#
# [pc40583@ppserver trimmed]$ ls -l ERR019315_*
# -rw-r--r-- 1 pc40583 ics       8749 Oct  2 10:46 ERR019315_1P.fastq.gz
# -rw-r--r-- 1 pc40583 ics 1000791110 Oct  2 10:47 ERR019315_1U.fastq.gz
# -rw-r--r-- 1 pc40583 ics       5089 Oct  2 10:46 ERR019315_2P.fastq.gz
# -rw-r--r-- 1 pc40583 ics        862 Oct  2 10:46 ERR019315_2U.fastq.gz

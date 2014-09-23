all: bloom_filter mapping

all_mito.fas: ../raw_data/references.txt
	python -c "import os, sys; sys.exit(os.system('cat %s > all_mito.fas' % ' '.join('../raw_data/%s.fasta' % f.strip() for f in open('../raw_data/references.txt'))))"

all_mito2.fas: all_mito.fas
	python -c "from Bio import SeqIO; SeqIO.write((r+r for r in SeqIO.parse('all_mito.fas', 'fasta')), 'all_mito2.fas', 'fasta')"

bloom_filter: all_mito2.fas
	../scripts/bloom_filter_all_mito.sh

mapping: all_mito2.fas bloom_filter
	../scripts/map_all_mito.sh

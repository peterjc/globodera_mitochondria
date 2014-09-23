all: bloom_filter mapping

all_mito2.fas: ../raw_data/references.txt
	python -c "import os, sys; sys.exit(os.system('cat %s > all_mito2.fas' % ' '.join('../raw_data/%s.fasta' % f.strip() for f in open('../raw_data/references.txt'))))"

bloom_filter: all_mito2.fas
	../scripts/bloom_filter_all_mito.sh

mapping: all_mito2.fas bloom_filter
	../scripts/map_all_mito.sh

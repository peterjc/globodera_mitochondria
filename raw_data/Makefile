all: download

download:
# Download files...
	wget -nc -c -i urls.txt
# "Making files read only...
	chmod a-w *.gb
	chmod a-w *.fasta
	chmod a-w *.fastq.gz
	echo "Download Done"

check:
# Check references.txt FASTA and GenBank files exist.
# Python snippet returns the number of missing files,
# thus zero is the desired non-error return value.
	python -c "import os, sys; sys.exit(sum(not os.path.isfile(f.strip()+'.gb') for f in open('references.txt')))"
	python -c "import os, sys; sys.exit(sum(not os.path.isfile(f.strip()+'.fasta') for f in open('references.txt')))"
# Check MD5 matches...
	md5sum -c md5sum.txt
	echo "Check done"

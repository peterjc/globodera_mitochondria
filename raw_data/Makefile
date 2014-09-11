all: ena

ena:
	wget -nc -c -i urls.txt
	chmod a-w *.gb
	chmod a-w *.fasta
	chmod a-w *.fastq.gz

check:
	md5sum -c md5sum.txt

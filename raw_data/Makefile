all: download

download:
# Download files...
	wget -nc -c -i urls.txt
	wget -O ERP000297.tsv "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERP000297&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,library_layout,nominal_length,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,submitted_bytes,submitted_md5,submitted_ftp,submitted_format,sra_bytes,sra_md5,sra_ftp,col_tax_id,col_scientific_name&download=text"
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
	python ../scripts/update_ENA_checksums.py ERP000297.tsv md5sum.txt
	md5sum -c md5sum.txt
	echo "Check done"

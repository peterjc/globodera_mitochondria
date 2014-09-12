#!/usr/bin/env python
import os
import sys
import urllib

project = "ERP000297"

#This old URL no longer works
#submissions_url = "http://www.ebi.ac.uk/ena/data/view/reports/sra/submitted_files/internal/%s" % project
#This URL picks all columns except the Aspera and Galaxy download links.
fields = "study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,library_layout,nominal_length,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,submitted_bytes,submitted_md5,submitted_ftp,submitted_format,sra_bytes,sra_md5,sra_ftp,col_tax_id,col_scientific_name".split(",")
project_url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run&fields=%s&download=text" % (project, ",".join(fields))
project_file = "../raw_data/%s.tsv" % project


def download_in_one(url, filename):
    print "Fetching %s" % url
    n = urllib.urlopen(url)
    data = n.read()
    n.close()

    h = open(filename, "w")
    h.write(data)
    h.close()
    print "Saved as %s" % filename

if not os.path.isfile(project_file):
    download_in_one(project_url, project_file)

def process_fastq(project, project_filename):
    h = open(project_filename)
    line = h.readline()
    #These were the columns in the old tabular output...
    #assert line == 'Study\tSample\tExperiment\tRun\tOrganism\tInstrument Platform\tInstrument Model\tLibrary Name\tLibrary Layout\tLibrary Source\tLibrary Selection\tRun Read Count\tRun Base Count\tFile Name\tFile Size\tmd5\tFtp\n', repr(line)
    assert line == "\t".join(fields) + "\n", repr(line)
    PRJ = fields.index("secondary_study_accession")
    URL = fields.index("fastq_ftp")
    for line in h:
        parts = line.rstrip("\n").split("\t")
        assert parts[PRJ] == project
        for url in parts[URL].split(";"):
            assert url.startswith("ftp.sra.ebi.ac.uk/vol1/fastq/"), url
            filename = os.path.join(os.path.dirname(project_filename), os.path.basename(url))
            if os.path.isfile(filename):
                print "Already have %s" % filename
                continue
            if not filename.endswith(".fastq.gz"):
                print "Skipping %s" % filename
                continue
            print "Downloading %s --> %s" % (url, filename)
            continue
            #Download file...
            rc = os.system("wget -O %s %s" % (filename, url))
            assert not rc, rc
            #Now check the md5...
            print filename
    h.close()

process_fastq(project, project_file)

#!/usr/bin/env python
import os
import sys

assert len(sys.argv) == 3
ena_table_filename = sys.argv[1]
checksum_filename = sys.argv[2]

def load_ENA_fastq_checksums(filename):
    handle = open(filename)
    headers = handle.readline().rstrip("\n").split("\t")
    fastq_url = headers.index("fastq_ftp")
    fastq_md5 = headers.index("fastq_md5")
    results = dict()
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        urls = fields[fastq_url].split(";")
        md5s = fields[fastq_md5].split(";")
        assert len(urls) == len(md5s)
        for url, md5 in zip(urls, md5s):
            filename = os.path.basename(url)
            results[filename] = md5
    handle.close()
    return results


def update_checksum_list(filename, md5dict):
    if os.path.isfile(filename):
        old = dict()
        with open(filename) as handle:
            for line in handle:
                if not line.strip():
                    continue
                md5, f = line.rstrip().split()
                old[f] = md5
        for k, v in md5dict.items():
            if k in old:
                assert v == old[k], "Mismatch for %s" % k
            else:
                old[k] = v
    else:
        old = md5dict
    with open(filename, "w") as handle:
        for f, md5 in sorted(old.items()):
            handle.write("%s  %s\n" % (md5, f))
    print("Now have %i MD5 checksums" % len(old))

md5dict = load_ENA_fastq_checksums(ena_table_filename)
print("Loaded checksums for %i ENA files" % len(md5dict))
update_checksum_list(checksum_filename, md5dict)

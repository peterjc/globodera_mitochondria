#!/usr/bin/env python
"""Find k-mers common to the known mtDNA circles.

This is a recreation of an old one-off script.

We find the following 53-mer sequences is present exactly in nice of
references (AJ249395, DQ631911, DQ631912, DQ631913, EF462976, EF462977,
EF462978, EF462980, EF462981). Missing in DQ631914, EF193005, EF462979.

'TTTGGTGGGTATTGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGA' or
'TCCATCAATCTGAACCTATATCATAGTCCCTCTGAGACTCAATACCCACCAAA' (rev. comp.)

This falls in the non-coding '222' region which shares similarity
between the circles. While it is not universal to the known circles,
it may make an interesting seed region for iterative assembly, or
be useful for wet-lab fishing for novel circles (e.g. with a PCR
primer to as a hybridisation target).
"""

from Bio.Seq import reverse_complement
from Bio import SeqIO

with open("../raw_data/references.txt") as handle:
    references = [x.strip() for x in handle]
fasta_files = ["../raw_data/%s.fasta" % acc for acc in references]
#print references

def get_kmers_from_circular_fasta(filename, kmer_length):
    """Takes a single-record FASTA which is treated as circular."""
    seq = str(SeqIO.read(filename, "fasta").seq).upper()
    length = len(seq)
    if length < kmer_length:
        raise StopIteration
    seq += seq[:kmer_length]
    for i in range(length):
        #kmer = seq[i:i+kmer_length]
        #assert len(kmer) == kmer_length
        #yield kmer
        yield seq[length:length+kmer_length]

def count_kmers_by_acc(circular_refs, kmer_length):
    if isinstance(circular_refs, str):
        raise TypeError("Want list of reference accessions")
    counts = dict()
    for acc in circular_refs:
        fasta = "../raw_data/%s.fasta" % acc
        for kmer in get_kmers_from_circular_fasta(fasta, kmer_length):
            if kmer in counts:
                counts[kmer].add(acc)
            else:
                rc = reverse_complement(kmer)
                if rc in counts:
                    counts[rc].add(acc)
                else:
                    counts[kmer] = set([acc])
    return counts

kmer_len = 25
min_count = 9
while True:
    print("-" * 50)
    print("k-mer length %i" % kmer_len)
    counts = count_kmers_by_acc(references, kmer_len)
    print("Found %i unique %i-mers in the %i references" % (len(counts), kmer_len, len(references)))
    match = False
    for target in range(len(references), min_count - 1, -1):
        match = False
        print("Present in %i references:" % target)
        for kmer, accs in sorted(counts.items()):
            if len(accs) == target:
                match = True
                print("%s in %s" % (kmer, ";".join(sorted(accs))))
        if match:
            break
    if match:
        kmer_len += 1
    else:
        break
print("Done")

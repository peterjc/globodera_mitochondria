Analysis of the mitochondria of *Globodera* nematodes
=====================================================

This repository is a work in progress (with the goal of forming part of a
scientific paper) using publicly released sequence data, primarily Illumina
reads for species *Globodera pallida* and *G. rostochiensis* from the BBSRC
funded *Globodera pallida* project (Cotton et al. 2014):

- Sequencing project: http://www.sanger.ac.uk/resources/downloads/helminths/globodera-pallida.html
- Sanger and 454 reads: ftp://ftp.sanger.ac.uk/pub/project/pathogens/Globodera/pallida/
- Illumina reads: http://www.ebi.ac.uk/ena/data/view/PRJEB2242

The mitochondria of these potato cyst nematodes are unusual for animals in
that rather than a single circular molecule of mitochondrial DNA (mtDNA),
they have multiple separate smaller circular molecules which together seem
to contain the expected gene complement although with duplication and
pseudogenes. These circles are known to vary between populations and can be
used as markers.


Published *Globodera* mtDNA circles
===================================

TODO - Add list here plus associated references.


Mapping & Coverage
==================

Few if any mainstream NGS mappers understand circular references, therefore
we are mapping to 'doubled' references for the mitochondrial genomes. Script
https://github.com/peterjc/picobio/blob/master/sambam/sam_circular_coverage.py
can then be used for calculating coverage modulo the real circle length.
[Currently this is a simplistic coverage calculation which does not skip
based jumped by a CIGAR D/N operator]. This script gives each read weight
1/n if the read was mapped to n locations.

Furthermore, we can 'undouble' the mapping for a more normal SAM/BAM files,
removing obvious duplicates if multiple mappings where recorded, albeit one
where reads can map off the right hand edge, using script
https://github.com/peterjc/picobio/blob/master/sambam/sam_undouble_circles.py

I have been trying mappers requesting all read mappings, and currently do
the coverage calculation on a per-circle basis. The coverage weighting
takes care of the need to de-duplicate so the results are near identical
with and without running ``sam_undouble_circles.py``.

However, requesting all mappings may introduce some artefacts near the origin.
e.g. Given a 100bp read which maps over the origin, you might get a partial
mapping reported at POS 0, CIGAR 90M, while the full mapping is at circle
length - 10, CIGAR 100M. These are not identical so deduplication will not
remove them - however the weighted coverage approach should handle this.

The key consideration with multiple mapping is placing reads which can be
mapped to multiple circles... perhaps exploiting partner reads here to flag
potential reads from similar but unknown mitochondrial circles?


TODO
====

Replace checked in ``scripts/blooming_reads.py`` etc with those already under
https://github.com/peterjc/picobio/tree/master/blooming_reads


References
==========

Cotton et al. (2004) The genome and life-stage specific transcriptomes of
*Globodera pallida* elucidate key aspects of plant parasitism by a cyst
nematode. *Genome Biology*. http://dx.doi.org/10.1186/gb-2014-15-3-r43

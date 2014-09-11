Raw data (read only once downloaded)
====================================

This folder holds raw data (e.g. Illumina reads from ENA) with scripts
to download and verify them (with a checksum), and thereafter they are
treated as read only for downstream analysis.

Note that due to its size, little (if any) of this raw data will be
checked into git - just the data fetching scripts and checksums will be.

To fetch the data (repeating this should be safe):

.. sourcecode:: console

    $ make

If you want to store the files elsewhere (or have previously downloaded
them) then using a manually created symlink should be fine.

To verify the data, which is done using md5 checksums (slow):

.. sourcecode:: console

    $ make check

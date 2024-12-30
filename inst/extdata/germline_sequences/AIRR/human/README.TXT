This is the extdata/germline_sequences/AIRR/human/ folder in the
igblastr package.

This folder contains the germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for Human.
Note that these FASTA files were obtained programmatically by running the
following code in the folder on Dec 29, 2024:

    igblastr:::download_human_germline_sequences_from_OGRDB()
    igblastr:::download_human_germline_sequences_from_OGRDB(extended=TRUE,
                                                            destdir="extended")

In addition to the FASTA files, this folder has a 'version' file that
contains the date of the download in YYYYMM format.


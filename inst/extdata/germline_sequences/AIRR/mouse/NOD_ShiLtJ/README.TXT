This is the extdata/germline_sequences/AIRR/mouse/NOD_ShiLtJ/ folder in the
igblastr package.

This folder contains the germline sequences (in FASTA format) downloaded
from AIRR-community/OGRDB for mouse strain NOD/ShiLtJ.
Note that these FASTA files were obtained programmatically by running the
following code in the folder on Dec 29, 2024:

    set_names <- paste("NOD/ShiLtJ", c("IGH", "IGKV", "IGLV"))
    igblastr:::download_mouse_germline_sequences_from_OGRDB(set_names)

In addition to the FASTA files, this folder has a 'version' file that
contains the date of the download in YYYYMM format.


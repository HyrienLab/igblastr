This is the extdata/germline_sequences/AIRR/mouse/LEWES_EiJ/ folder in the
igblastr package.

This folder contains the germline sequences (in FASTA format) downloaded
from AIRR-community/OGRDB for mouse strain LEWES/EiJ.
Note that these FASTA files were obtained programmatically by running the
following code in the folder on Dec 29, 2024:

    set_names <- c("LEWES/EiJ IGH", "LEWES/EiJ IGKV", "LEWES/EiJ IGLV")
    igblastr:::download_mouse_germline_sequences_from_OGRDB(set_names)


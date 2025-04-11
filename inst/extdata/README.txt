This is the extdata/ folder in the igblastr package.

This folder contains various data files used in the man page examples and/or
unit tests of the igblastr package.

Content
=======

- README.txt: This file.

- 1279067_1_Paired_sequences.fasta: FASTA file containing 8437 pairs of
  human antibody sequences (16874 individual sequences) retrieved from OAS (the
  Observed Antibody Space database). The file was obtained programmatically
  by running the following code in this folder on March 26, 2025:

    library(igblastr)
    download_paired_OAS_units("Jaffe_2022", "1279067_1_Paired_All.csv.gz")
    df <- read_OAS_csv("1279067_1_Paired_All.csv.gz")
    sequences <- extract_sequences_from_paired_OAS_df(df, add.prefix=TRUE)
    writeXStringSet(sequences, "1279067_1_Paired_sequences.fasta")

- 1279067_1_Paired_All.json: JSON file containing the metadata associated
  with the sequences in 1279067_1_Paired_sequences.fasta.
  This file was obtained by downloading 1279067_1_Paired_All.json directly
  from https://opig.stats.ox.ac.uk/webapps/ngsdb/paired/Jaffe_2022/json/

- catnap_bnabs.fasta: FASTA file containing 1000 heavy- and light-chains
  associated with bnAbs downloaded from the CATNAP database.

- constant_regions/: See README.txt in constant_regions/ subfolder.

- germline_sequences/: See README.txt in germline_sequences/ subfolder.

- ncbi_igblast_data_files/: See README.txt in ncbi_igblast_data_files/
  subfolder.


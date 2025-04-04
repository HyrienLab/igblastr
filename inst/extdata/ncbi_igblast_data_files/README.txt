This is the extdata/ncbi_igblast_data_files/ folder in the igblastr package.

This folder contains subfolders internal_data/ and optional_file/ extracted
from the source tarball ncbi-igblast-1.22.0-src.tar.gz provided at
https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/

More precisely, internal_data/ and optional_file/ were
obtained by:
  1. Downloading ncbi-igblast-1.22.0-src.tar.gz from
       https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/
  2. Running:
       tar zxf ncbi-igblast-1.22.0-src.tar.gz
       mv ncbi-igblast-1.22.0-src/c++/src/app/igblast/internal_data .
       mv ncbi-igblast-1.22.0-src/c++/src/app/igblast/optional_file .

-----------------------------------------------------------------------------

Why are internal_data/ and optional_file/ included in the igblastr package?

On a Mac Silicon machine, the install_igblast() function in igblastr uses
disk image ncbi-igblast-1.22.0+.dmg available at NCBI to install IgBLAST.
However, for some mysterious reason, this disk image does NOT include
the internal_data/ or optional_file/ folders so install_igblast() fixes
this by copying them to the _root directory_ of the IgBLAST installation.


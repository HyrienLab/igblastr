**igblastr** is an R package that provides functions to conveniently install
and use a local IgBLAST installation from R.

IgBLAST is described at https://pubmed.ncbi.nlm.nih.gov/23671333/

Online IgBLAST: https://www.ncbi.nlm.nih.gov/igblast/

Note that the package is still a WORK-IN-PROGRESS! Please use
https://github.com/hpages/igblastr/issues to report bugs, provide
feedback, request features, etc...


## Quick start


### 1. Install and load igblastr

#### Install Bioconductor dependencies

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install(c('Biostrings','S4Vectors','GenomeInfoDb'))

#### Install igblastr

    if (!require(remotes, quietly=TRUE))
        install.packages("remotes")

    library(remotes)
    install_github("hpages/igblastr")

#### Load igblastr

    library(igblastr)


### 2. Set up IgBLAST

#### Install a pre-compiled IgBLAST

    install_igblast()

This will download the arch-specific pre-compiled IgBLAST from
https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/
and install it in a persistent location. See `?install_igblast`.

#### or use an existing installation of IgBLAST

    Sys.setenv(IGBLAST_ROOT="path/to/igblast/root")

Note: Environment variable `IGBLAST_ROOT` should preferrably be set
in a more persistent manner outside R e.g. by defining it in the
user's `.profile` (on Unix/Mac). See `?IGBLAST_ROOT`.

#### Check IgBLAST

    igblast_info()

This will make sure that the installed IgBLAST is functional. It will
also display basic information about it (e.g. `igblastn` version).


### 3. Install germline dbs and select db to use

#### Install germline dbs

At least one germline db must be installed.

Several specialized functions will be provided for that e.g.:

    install_NCBI_germline_db()  # not ready yet!
    install_IMGT_germline_db()  # see '?install_IMGT_germline_db'
    install_AIRR_germline_db()  # not ready yet!
    etc...

Note:

- `install_NCBI_germline_db()` will manage installation of the dbs available
  at https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/.
  These are:
  - `ncbi_human_c_genes`
  - `mouse_gl_VDJ`
  - `rhesus_monkey_VJ`
  - `airr_c_human`
  - `airr_c_mouse`

#### Select the db to use with igblastn()

    list_germline_dbs()
    use_germline_db(<db-name>)

See `?use_germline_db`.


### 4. Use igblastn()

See `?igblastn`.


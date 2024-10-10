WORK IN PROGRESS!


## Tentative workflow


4 easy steps



### 1. Install and load igblastr


#### Install

    library(remotes)
    install_github("hpages/igblastr")


#### Load

    library(igblastr)



### 2. Set up IgBlast


#### Use an existing installation of IgBlast

    Sys.setenv(IGBLAST_ROOT, "path/to/igblast/root")

Note: Environment variable IGBLAST\_ROOT should preferrably be set
in a more persistent fashion outside R e.g. by defining it in the
user's `.profile` (on Unix/Mac).


#### or install a pre-compiled IgBlast binary

    install_igblast()

This will download the arch-specific tarball from
https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/
and extract it in a persistent user-writable hidden location (e.g.
under `R_user_dir("igblastr", "cache")`).

Notes:
- No tarball for arm64 Mac (Mac Silicon), only a `dmg`.
- This installation takes precedence over the installation indicated
  by IGBLAST\_ROOT.


#### Check IgBlast

    igblast_info()

This makes sure that the installed IgBlast is functional. Also displays
some info about it (e.g. `igblastn` and `igblastp` versions).



### 3. Install germline dbs and select db to use


#### Install germline dbs

At least one germline db must be installed.

Several specialized functions will be provided for that:

    install_NCBI_germline_db()
    install_VQUEST_germline_db()
    install_OGRDB_germline_db()
    install_AIRR_germline_db()
    etc...

`install_NCBI_germline_db()` will manage installation of the dbs available
at https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/.

These are:
- `ncbi_human_c_genes`
- `mouse_gl_VDJ`
- `rhesus_monkey_VJ`
- `airr_c_human`
- `airr_c_mouse`

`install_VQUEST_germline_db()` will manage installation of the dbs available
at https://www.imgt.org/download/V-QUEST/. It will require Perl!
The following naming scheme will be used for the names of the installed dbs:
`VQUEST_<version>_<organism>_<IG|TR|full>`.


#### List the installed germline dbs

    list_germline_dbs()

This will produce a `data.frame` with 1 row per installed db:
- 1st col: name of the db
- 2nd col: organism in binomial (latin) species form e.g. `Homo sapiens`


#### Select the db to be used with igblastn() or igblastp()

    use_germline_db(<db-name>)



### 4. Use igblastn() or igblastp()


Coming soon...


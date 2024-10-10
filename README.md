WORK IN PROGRESS!


## Tentative workflow


    library(igblastr)



### Set up IgBlast


#### Use an existing installation of IgBlast

    Sys.setenv(IGBLAST_ROOT, "path/to/igblast/root")
    
    # Note: Environment variable IGBLAST_ROOT should preferrably be set once
    #       for all outside R e.g. by defining it in the user's .profile
    #       on Unix/Mac.


#### or install a pre-compiled IgBlast binary

    install_igblast()
    
    # --> Downloads the arch-specific tarball from
    #     https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/
    #     and extracts under R_user_dir("igblastr", "cache")
    # Note: No tarball for arm64 Mac (Mac Silicon), only a dmg.
    # Note: This installation takes precedence over the installation indicated
    #       by IGBLAST_ROOT.


#### Check IgBlast

    igblast_info()
    
    # --> Makes sure that the installed IgBlast is functional and displays
    #     some info about it (e.g. igblastn and igblastp versions).



### Select germline db

#### List the installed germline dbs

    list_germline_dbs()
    
    # --> Produces a data.frame with 1 row per db.
    #     1st col is the name of the db, 2nd col is the organism (binomial
    #     species e.g. Homo sapiens).
    #
    #     Pre-installed dbs (from NCBI):
    #       ncbi_human_c_genes
    #       mouse_gl_VDJ
    #       rhesus_monkey_VJ
    #       airr_c_human
    #       airr_c_mouse
    #     See below for more dbs.


#### Install additional germline dbs if necessary

Several specialized functions will be provides e.g.

    install_VQUEST_germline_db()  # requires Perl
    
    # --> Will install dbs named as:
    #     VQUEST_<version>_<organism>_<IG|TR|full>

or

    install_OGRDB_germline_db()
    install_AIRR_germline_db()
    etc...


#### Select the db to be used with igblastn() or igblastp()

    use_germline_db(<db-name>)



### Use igblastn()
    
    coming soon...


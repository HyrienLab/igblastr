### =========================================================================
### reset_germline_dbs()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .install_missing_builtin_OGRDB_human_germline_dbs()
###

.form_builtin_OGRDB_human_germline_db_name <- function(fasta_store)
{
    stopifnot(isSingleNonWhiteString(fasta_store))
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- basename(dirname(fasta_store))
    db_name <- sprintf("_OGRDB.human.%s.%s", loci_in1string, version)
    flavor <- basename(fasta_store)
    stopifnot(flavor %in% c("ref", "src"))
    if (flavor == "src")
        db_name <- paste0(db_name, ".src")
    db_name
}

### Install db only if missing.
.install_builtin_OGRDB_human_germline_db <-
    function(install_dir, fasta_store, verbose=FALSE)
{
    db_name <- .form_builtin_OGRDB_human_germline_db_name(fasta_store)
    install_germline_db(install_dir, db_name, fasta_store, IG_LOCI,
                        gapped=TRUE, with.intdata=TRUE,
                        if.exists="no-op", verbose=verbose)
    db_path <- file.path(install_dir, db_name)
    add_computed_auxdata_to_germline_db(db_path, verbose=verbose)
}

.install_missing_builtin_OGRDB_human_germline_dbs <-
    function(install_dir, human_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(human_dir), dir.exists(human_dir),
              isTRUEorFALSE(verbose))
    subdirs <- list.dirs(human_dir, full.names=FALSE, recursive=FALSE)
    ## Starting with igblastr 1.0.21/1.1.21, versions 202309 and 202401
    ## are excluded from the list of built-in OGRDB germline dbs for human.
    ## The reason for this is that we found out that the internal
    ## data included in OGRDB germline set IGLambda_VJ versions 1 & 2
    ## is inconsistent. See validate_human_intdata() in
    ## extdata/germline_sets/OGRDB/human/download_human_germline_sequences.R
    ## for more information.
    exclude_list <- c("202309", "202401")
    subdirs <- setdiff(subdirs, c("diffs", exclude_list))
    fasta_stores <- file.path(human_dir, rep(subdirs, each=2L), c("ref", "src"))
    for (fasta_store in fasta_stores)
        .install_builtin_OGRDB_human_germline_db(install_dir, fasta_store,
                                                 verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .install_missing_builtin_OGRDB_mouse_germline_dbs()
###

.form_builtin_OGRDB_mouse_germline_db_name <- function(fasta_store)
{
    stopifnot(isSingleNonWhiteString(fasta_store))
    strain <- basename(fasta_store)
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- read_version_file(fasta_store)
    sprintf("_OGRDB.mouse.%s.%s.%s", strain, loci_in1string, version)
}

### Install db only if missing.
.install_builtin_OGRDB_mouse_germline_db <-
    function(install_dir, fasta_store, verbose=FALSE)
{
    db_name <- .form_builtin_OGRDB_mouse_germline_db_name(fasta_store)
    install_germline_db(install_dir, db_name, fasta_store, IG_LOCI,
                        gapped=TRUE, with.intdata=TRUE,
                        if.exists="no-op", verbose=verbose)
    ## The computed data for mouse strains NOD_ShiLtJ and PWD_PhJ version
    ## 202501 is incomplete (i.e. contains NAs) so there's no need to call
    ## add_computed_auxdata_to_germline_db(). If we do so, then the auxdata
    ## will be computed and discarded each time list_germline_dbs() is called!
    #db_path <- file.path(install_dir, db_name)
    #add_computed_auxdata_to_germline_db(path, verbose=verbose)
}

.install_missing_builtin_OGRDB_mouse_germline_dbs <-
    function(install_dir, mouse_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(mouse_dir), dir.exists(mouse_dir),
              isTRUEorFALSE(verbose))
    fasta_stores <- list.dirs(mouse_dir, recursive=FALSE)
    ## Starting with igblastr 1.0.21/1.1.21, mouse strains CAST/EiJ, LEWES/EiJ,
    ## and MSM/MsJ are excluded from the list of built-in OGRDB germline dbs
    ## for mouse. The reason for this is that we found out that the internal
    ## data included in OGRDB germline sets "CAST/EiJ IGH" version 1,
    ## "CAST/EiJ IGLV" version 1, "LEWES/EiJ IGH" version 1, and
    ## "MSM/MsJ IGLV" version 1 is inconsistent. See validate_mouse_intdata()
    ## in extdata/germline_sets/OGRDB/mouse/download_mouse_germline_sequences.R
    ## for more information.
    exclude_list <- c("CAST_EiJ", "LEWES_EiJ", "MSM_MsJ")
    fasta_stores <- fasta_stores[!(basename(fasta_stores) %in% exclude_list)]
    for (fasta_store in fasta_stores)
        .install_builtin_OGRDB_mouse_germline_db(install_dir, fasta_store,
                                                 verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .install_missing_builtin_OGRDB_rhesus_monkey_germline_dbs()
###

.form_builtin_OGRDB_rhesus_monkey_germline_db_name <- function(fasta_store)
{
    stopifnot(isSingleNonWhiteString(fasta_store))
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- basename(fasta_store)
    sprintf("_OGRDB.rhesus_monkey.%s.%s", loci_in1string, version)
}

### Install db only if missing.
.install_builtin_OGRDB_rhesus_monkey_germline_db <-
    function(install_dir, fasta_store, verbose=FALSE)
{
    db_name <- .form_builtin_OGRDB_rhesus_monkey_germline_db_name(fasta_store)
    install_germline_db(install_dir, db_name, fasta_store, IG_LOCI,
                        gapped=TRUE, with.intdata=TRUE,
                        if.exists="no-op", verbose=verbose)
    db_path <- file.path(install_dir, db_name)
    add_computed_auxdata_to_germline_db(db_path, verbose=verbose)
}

.install_missing_builtin_OGRDB_rhesus_monkey_germline_dbs <-
    function(install_dir, rhesus_monkey_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(rhesus_monkey_dir),
              dir.exists(rhesus_monkey_dir), isTRUEorFALSE(verbose))
    subdirs <- list.dirs(rhesus_monkey_dir, full.names=FALSE, recursive=FALSE)
    ## Starting with igblastr 1.0.21/1.1.21, we replaced version 202601 with
    ## version 202602.
    exclude_list <- "202601"
    subdirs <- setdiff(subdirs, c("diffs", exclude_list))
    fasta_stores <- file.path(rhesus_monkey_dir, subdirs)
    for (fasta_store in fasta_stores)
        .install_builtin_OGRDB_rhesus_monkey_germline_db(install_dir,
                                                         fasta_store,
                                                         verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_missing_builtin_germline_dbs()
###

### Not exported!
### 'install_dir' must exist.
install_missing_builtin_germline_dbs <- function(install_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(install_dir), dir.exists(install_dir),
              isTRUEorFALSE(verbose))

    OGRDB_germline_seq_dir <- system.file(package="igblastr",
                                 "extdata", "germline_sets", "OGRDB",
                                 mustWork=TRUE)

    human_dir <- file.path(OGRDB_germline_seq_dir, "human")
    .install_missing_builtin_OGRDB_human_germline_dbs(install_dir, human_dir,
                                                      verbose=verbose)

    mouse_dir <- file.path(OGRDB_germline_seq_dir, "mouse")
    .install_missing_builtin_OGRDB_mouse_germline_dbs(install_dir, mouse_dir,
                                                      verbose=verbose)

    rhesus_monkey_dir <- file.path(OGRDB_germline_seq_dir, "rhesus_monkey")
    .install_missing_builtin_OGRDB_rhesus_monkey_germline_dbs(install_dir,
                                                              rhesus_monkey_dir,
                                                              verbose=verbose)

    ## Any other built-in germline dbs to install?
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reset_germline_dbs()
###

### Has its own tests in tests/testthat/test-reset_germline_dbs.R!
### If 'install_dir' exists, it gets destroyed and replaced with a freshly
### populated directory. In other words, its final content doesn't depend
### on whether it already exists or not. Note that the destroy/replace
### operation is atomic.
.install_all_builtin_germline_dbs <- function(install_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(install_dir))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'install_dir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_install_dir <- tempfile("builtin_germline_dbs_")
    dir.create(tmp_install_dir)
    on.exit(nuke_file(tmp_install_dir))

    install_missing_builtin_germline_dbs(tmp_install_dir, verbose=verbose)

    ## Everything went fine so we can rename 'tmp_install_dir' to 'install_dir'.
    rename_file(tmp_install_dir, install_dir, replace=TRUE)
}

### Will nuke any user-installed db!
reset_germline_dbs <- function(verbose=FALSE)
{
    set_db_in_use("germline", "")  # cancel current selection
    germline_dbs_home <- igblastr_cache(GERMLINE_DBS)
    .install_all_builtin_germline_dbs(germline_dbs_home, verbose=verbose)
}


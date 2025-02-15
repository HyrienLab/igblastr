### =========================================================================
### create_c_region_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.list_C_fasta_files <- function(dirpath)
{
    stopifnot(isSingleNonWhiteString(dirpath))
    pattern <- paste0("^IG[HKL]C\\.fasta$")
    files <- list.files(dirpath, pattern=pattern)
    if (length(files) == 0L)
        stop(wmsg("Anomaly: no C-region sequence files found in ", dirpath))
    files
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_c_region_db()
###

.stop_on_existing_c_region_db <- function(destdir)
{
    db_name <- basename(destdir)
    msg1 <- c("C-region db ", db_name, " is already installed.")
    msg2 <- c("Use list_c_region_dbs() to list the C-region databases ",
              "already installed on your machine (see '?list_c_region_dbs').")
    msg3 <- c("Use 'force=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### Perl required!
###
### Creates a C-region db (constant regions) from a collection of FASTA
### files (typically obtained from IMGT) for a given organism.
### Note that 'destdir' will typically be the path to a subdir of the
### C_REGION_DBS cache compartment (see R/cache-utils.R for details about
### igblastr's cache organization). This subdir or any of its parent
### directories don't need to exist yet.
create_c_region_db <- function(fastadir, destdir, force=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (dir.exists(destdir) && !force)
        .stop_on_existing_c_region_db(destdir)
    edit_fasta_script <- get_edit_imgt_file_Perl_script()
    fasta_files <- .list_C_fasta_files(fastadir)
    fasta_files <- file.path(fastadir, fasta_files)

    ## We first create the db in a temporary folder, and, only if successful,
    ## replace 'destdir' with the temporary folder. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity and avoids
    ## loosing the content of the existing 'destdir' in case something goes
    ## wrong.
    tmp_destdir <- tempfile("c_region_db_")
    dir.create(tmp_destdir, recursive=TRUE)
    on.exit(nuke_file(tmp_destdir))
    create_region_db(fasta_files, tmp_destdir, region_type="C",
                     edit_fasta_script=edit_fasta_script)
    rename_file(tmp_destdir, destdir, replace=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### form_IMGT_c_region_db_name()
###

form_IMGT_c_region_db_name <- function(organism_path)
{
    fasta_files <- .list_C_fasta_files(organism_path)
    loci <- paste(sort(substr(fasta_files, 1L, 3L)), collapse="+")
    version <- read_version_file(organism_path)
    organism <- basename(organism_path)
    ## Prefix name with underscore because it's a builtin db.
    paste("_IMGT", organism, loci, version, sep=".")
}


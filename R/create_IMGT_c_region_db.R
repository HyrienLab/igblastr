### =========================================================================
### create_IMGT_c_region_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_IMGT_c_region_db()
###

.list_IMGT_C_fasta_files <- function(dirpath)
{
    pattern <- paste0("^IG[HKL]C\\.fasta$")
    files <- list.files(dirpath, pattern=pattern)
    if (length(files) == 0L)
        stop(wmsg("Anomaly: no C-region files found in ", dirpath))
    file.path(dirpath, files)
}

.stop_on_existing_IMGT_c_region_db <- function(destdir)
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
### Creates a C-region db (constant regions) from the FASTA files provided by
### IMGT for a given organism.
###
### 'organism_path' must be the path to one of the subfolders in
### 'system.file(package="igblastr", "extdata", "constant_regions", "IMGT")'.
###
### 'destdir' will typically be the path to a subdir of
###  <igblastr-cache>/c_region_dbs/
create_IMGT_c_region_db <- function(organism_path, destdir, force=FALSE)
{
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    edit_fasta_script <- get_edit_imgt_file_Perl_script()
    if (dir.exists(destdir) && !force)
        .stop_on_existing_IMGT_c_region_db(destdir)

    ## We first make the db in a temporary folder, and if successful, we
    ## replace 'destdir' with the temporary folder. This achieves atomicity
    ## and avoids loosing the content of the existing 'destdir' in case
    ## something goes wrong.
    tmpdestdir <- tempfile("c_region_db_")
    dir.create(tmpdestdir, recursive=TRUE)
    on.exit(unlink(tmpdestdir, recursive=TRUE, force=TRUE))
    fasta_files <- .list_IMGT_C_fasta_files(organism_path)
    create_IMGT_region_db(fasta_files, tmpdestdir, gene_segment="C",
                          edit_fasta_script=edit_fasta_script)
    replace_file(destdir, tmpdestdir)
}


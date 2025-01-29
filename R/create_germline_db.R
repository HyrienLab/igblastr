### =========================================================================
### create_germline_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### A thin wrapper around create_region_db().
.create_region_db2 <- function(fastadir, destdir, list_files_FUN,
                               edit_fasta_script, region_type=c("V", "D", "J"))
{
    region_type <- match.arg(region_type)
    fasta_files <- list_files_FUN(fastadir)
    create_region_db(fasta_files, destdir, region_type=region_type,
                     edit_fasta_script=edit_fasta_script)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .list_V_fasta_files()
### .list_D_fasta_files()
### .list_J_fasta_files()

.list_fasta_files <- function(dirpath, region_type=c("V", "D", "J"),
                              expected_files)
{
    region_type <- match.arg(region_type)
    pattern <- paste0(region_type, "\\.fasta$")
    files <- list.files(dirpath, pattern=pattern)
    if (length(files) == 0L)
        stop(wmsg("Anomaly: no ", region_type, " files found in ", dirpath))
    if (!setequal(files, expected_files)) {
        if (all(files %in% expected_files)) {
            warning("some ", region_type, " files are missing ",
                    "in ", dirpath, "/ compared to Homo_sapiens")
        } else {
            warning("set of ", region_type, " files in ", dirpath, "/ not ",
                    "the same as for Homo_sapiens")
        }
    }
    file.path(dirpath, files)
}

.list_V_fasta_files <- function(dirpath)
{
    EXPECTED_FILES <- paste0("IG", c("HV", "KV", "LV"), ".fasta")
    .list_fasta_files(dirpath, "V", EXPECTED_FILES)
}

.list_D_fasta_files <- function(dirpath)
{
    EXPECTED_FILES <- paste0("IGHD", ".fasta")
    .list_fasta_files(dirpath, "D", EXPECTED_FILES)
}

.list_J_fasta_files <- function(dirpath)
{
    EXPECTED_FILES <- paste0("IG", c("HJ", "KJ", "LJ"), ".fasta")
    .list_fasta_files(dirpath, "J", EXPECTED_FILES)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_germline_db()
###

.stop_on_existing_germline_db <- function(destdir)
{
    db_name <- basename(destdir)
    msg1 <- c("Germline db ", db_name, " is already installed.")
    msg2 <- c("Use list_germline_dbs() to list the germline databases ",
              "already installed on your machine (see '?list_germline_dbs').")
    msg3 <- c("Use 'force=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### Perl required!
###
### A "germline db" is made of three "region dbs": one V-, one D-, and one
### J-region db. Calls create_region_db() to create each "region db".
###
### 'destdir' will typically be the path to a subdir of:
###   <igblastr-cache>/germline_dbs/
### This subdir or any of its parent directories don't need to exist yet.
create_germline_db <- function(fastadir, destdir, force=FALSE)
{
    stopifnot(isSingleNonWhiteString(fastadir),
              isSingleNonWhiteString(destdir))
    if (!dir.exists(fastadir))
        stop(wmsg("Anomaly: directory ", fastadir, " not found"))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (dir.exists(destdir) && !force)
        .stop_on_existing_germline_db(destdir)
    edit_fasta_script <- get_edit_imgt_file_Perl_script()

    ## We first make the db in a temporary folder, and, if successful, we
    ## replace 'destdir' with the temporary folder. This achieves atomicity
    ## and avoids loosing the content of the existing 'destdir' in case
    ## something goes wrong.
    tmp_destdir <- tempfile("germline_db_")
    dir.create(tmp_destdir, recursive=TRUE)
    on.exit(nuke_file(tmp_destdir))
    .create_region_db2(fastadir, tmp_destdir, .list_V_fasta_files,
                       edit_fasta_script, region_type="V")
    .create_region_db2(fastadir, tmp_destdir, .list_D_fasta_files,
                       edit_fasta_script, region_type="D")
    .create_region_db2(fastadir, tmp_destdir, .list_J_fasta_files,
                       edit_fasta_script, region_type="J")
    rename_file(tmp_destdir, destdir, replace=TRUE)
}


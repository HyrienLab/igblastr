### =========================================================================
### create_c_region_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.LIST_C_REGION_DB_TIP <- c(
    "Use list_c_region_dbs() to list all the C-region dbs ",
    "currently installed in the cache (see '?list_c_region_dbs')."
)

stop_on_existing_c_region_db <- function(db_name)
{
    msg1 <- c("C-region db ", db_name, " is already installed ",
              "in igblastr's persistent cache.")
    msg3 <- c("Use 'overwrite=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(.LIST_C_REGION_DB_TIP), "\n  ", wmsg(msg3))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_loci_in_c_region_fasta_dir()
###

.list_C_fasta_files <- function(fasta_dir, loci_prefix)
{
    stopifnot(isSingleNonWhiteString(fasta_dir), dir.exists(fasta_dir),
              isSingleNonWhiteString(loci_prefix))
    pattern <- paste0("^", loci_prefix, ".C\\.fasta$")
    fasta_files <- list_fasta_files(fasta_dir, pattern, full.names=FALSE)
    stopifnot(length(fasta_files) != 0L)
    fasta_files
}

### Returns a character vector of loci in canonical order.
.get_loci_from_c_region_fasta_set <- function(fasta_files, loci_prefix)
{
    stopifnot(is.character(fasta_files),
              isSingleString(loci_prefix), loci_prefix %in% c("IG", "TR"))
    loci <- unique(sub("C\\.fasta$", "", fasta_files))
    sort_unique_loci(loci, loci_prefix=loci_prefix)
}

list_loci_in_c_region_fasta_dir <- function(fasta_dir, loci_prefix)
{
    fasta_files <- .list_C_fasta_files(fasta_dir, loci_prefix)
    .get_loci_from_c_region_fasta_set(fasta_files, loci_prefix)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .collect_C_fasta_files()
###

.collect_C_fasta_files <- function(fasta_dir, loci)
{
    if (!isSingleNonWhiteString(fasta_dir))
        stop(wmsg("'fasta_dir' must be a single (non-empty) string"))
    if (!dir.exists(fasta_dir))
        stop(wmsg("directory ", fasta_dir, " not found"))
    check_selected_loci(loci)
    wanted_files <- file.path(fasta_dir, paste0(loci, "C.fasta"))
    stopifnot(all(file.exists(wanted_files)))
    wanted_files
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_c_region_db()
###

### Creates a C-region db (constant regions) from a collection of FASTA
### files (typically obtained from IMGT) for a given organism.
### Note that 'destdir' will typically be the path to a subdir of the
### C_REGION_DBS cache compartment (see R/cache-utils.R for details about
### igblastr's cache organization). This subdir or any of its parent
### directories don't need to exist yet.
create_c_region_db <- function(fasta_dir, loci, destdir,
                               disambiguate.allele.names=FALSE,
                               overwrite=FALSE, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir))
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    if (verbose) {
        what <- c("CREATING ", basename(destdir))
        message("\n====== START ", wmsg(what), " ======\n")
    }

    if (dir.exists(destdir) && !overwrite)
        stop_on_existing_c_region_db(basename(destdir))

    fasta_files <- .collect_C_fasta_files(fasta_dir, loci)

    ## We first create the db in a temporary folder, and, only if successful,
    ## we replace 'destdir' with the temporary folder. Otherwise we destroy
    ## the temporary folder and raise an error. This achieves atomicity and
    ## avoids loosing the content of the existing 'destdir' in case something
    ## goes wrong.
    tmp_destdir <- tempfile("c_region_db_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))
    create_region_db(tmp_destdir, fasta_files, region_type="C",
                     disambiguate.allele.names=disambiguate.allele.names,
                     verbose=verbose)
    rename_file(tmp_destdir, destdir, replace=TRUE)

    if (verbose)
        message("====== DONE ", wmsg(what), " ======\n")
}


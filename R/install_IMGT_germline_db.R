### =========================================================================
### install_IMGT_germline_db()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_effective_loci()
###

### Returns the intersection between 'wanted_loci' and 'found_loci'
### but with lots of sanity checks, bells and whistles.
.get_effective_loci <- function(wanted_loci, found_loci)
{
    stopifnot(is.character(wanted_loci), is.character(found_loci))
    keep_idx <- which(wanted_loci %in% found_loci)
    if (length(keep_idx) == 0L) {
        what <- if (length(wanted_loci) == 1L) "locus" else "loci"
        in1string <- paste0(wanted_loci, collapse=", ")
        stop(wmsg("no FASTA files found for ", what, " ", in1string))
    }
    missing_loci <- wanted_loci[-keep_idx]
    ## Like 'intersect(found_loci, wanted_loci)' but the returned
    ## intersection is guaranteed to be ordered like in 'found_loci'.
    loci <- found_loci[found_loci %in% wanted_loci]
    if (length(missing_loci) != 0L) {
        what1 <- if (length(missing_loci) == 1L) "locus" else "loci"
        in1string1 <- paste0(missing_loci, collapse=", ")
        what2 <- if (length(loci) == 1L) "locus" else "loci"
        in1string2 <- paste0(loci, collapse=", ")
        warning(wmsg("No FASTA files found for ", what1, " ", in1string1, " ",
                     "--> Installing germline db for ", what2, " ", in1string2,
                     "."),
                immediate.=TRUE)
    }
    loci
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_germline_db()
###

.form_IMGT_germline_db_name <- function(fasta_store, loci)
{
    stopifnot(isSingleNonWhiteString(fasta_store), dir.exists(fasta_store))
    checkarg_loci(loci)
    organism_path <- dirname(fasta_store)
    organism <- basename(organism_path)
    refdir <- dirname(organism_path)
    stopifnot(basename(refdir) == VQUEST_REFERENCE_DIRECTORY)
    IMGT_store <- dirname(refdir)
    release <- basename(IMGT_store)
    sprintf("IMGT-%s.%s.%s", release, organism, paste(loci, collapse="+"))
}

install_IMGT_germline_db <- function(release, organism="Homo sapiens",
                                     tcr.db=FALSE, loci="auto",
                                     without.intdata=FALSE,
                                     without.auxdata=FALSE,
                                     overwrite=FALSE, verbose=FALSE, ...)
{
    ## Check arguments.
    organism <- normalize_IMGT_organism(organism)
    loci <- normalize_user_supplied_loci(loci, tcr.db=tcr.db,
                                         stop.if.missing.regions=TRUE)
    loci_prefix <- extract_loci_prefix(loci)
    if (!isTRUEorFALSE(without.intdata))
        stop(wmsg("'without.intdata' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(without.auxdata))
        stop(wmsg("'without.auxdata' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## Download IMGT/V-QUEST release to local store if it's not already there.
    IMGT_store <- download_IMGT_release_to_IMGT_store(release, ...)

    ## Get path to local FASTA store for IMGT germline sequences.
    fasta_store <- find_organism_fasta_store_in_IMGT_store(IMGT_store,
                                                           organism,
                                                           loci_prefix)

    ## Keep loci for which IMGT actually provides FASTA files.
    found_loci <- list_loci_in_germline_fasta_dir(fasta_store, loci_prefix)
    loci <- .get_effective_loci(loci, found_loci)

    ## Compute 'db_name'.
    db_name <- .form_IMGT_germline_db_name(fasta_store, loci)

    ## Create and install germline db.
    install_dir <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    with.auxdata <- !(without.auxdata || loci_prefix == "TR")
    if.exists <- if (overwrite) "overwrite" else "error"
    install_germline_db(install_dir, db_name, fasta_store, loci,
                        gapped=TRUE, with.intdata=!without.intdata,
                        with.auxdata=with.auxdata,
                        if.exists=if.exists, verbose=verbose,
                        cheer.if.success=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_c_region_db()
###

.form_IMGT_c_region_db_name <- function(organism, loci, version)
{
    stopifnot(isSingleNonWhiteString(organism))
    sprintf("IMGT.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

install_IMGT_c_region_db <- function(organism, loci,
                                     disambiguate.allele.names=FALSE,
                                     overwrite=FALSE, verbose=FALSE)
{
    organism <- find_organism_shortname(normalize_IMGT_organism(organism))
    loci <- normalize_user_supplied_loci(loci)
    loci_prefix <- extract_loci_prefix(loci)
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## Get path to local FASTA store for IMGT C-region sequences.
    fasta_store <- path_to_IMGT_c_region_fasta_store(organism, loci_prefix)
    if (!dir.exists(fasta_store))
        stop(wmsg("we're not aware of any ", loci_prefix, " C-region ",
                  "sequences available at IMGT for ", organism, ", sorry"))

    ## Keep loci for which IMGT actually provides FASTA files.
    found_loci <- list_loci_in_c_region_fasta_dir(fasta_store, loci_prefix)
    loci <- .get_effective_loci(loci, found_loci)

    ## Compute 'db_name'.
    version <- read_version_file(fasta_store)
    db_name <- .form_IMGT_c_region_db_name(organism, loci, version)

    ## Create IMGT C-region db.
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    db_path <- file.path(c_region_dbs_home, db_name)
    create_c_region_db(fasta_store, loci, db_path,
                       disambiguate.allele.names=disambiguate.allele.names,
                       overwrite=overwrite, verbose=verbose)

    ## Success!
    message("C-region db ", db_name, " successfully installed.")
    message("Call use_c_region_db(\"", db_name, "\") to select it as the")
    message("C-region db to use with igblastn().")

    invisible(db_name)
}


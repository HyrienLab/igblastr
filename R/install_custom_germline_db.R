### =========================================================================
### install_custom_germline_db()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_germline_db()
###

.LIST_GERMLINE_DB_TIP <- c(
    "Use list_germline_dbs() to list all the germline dbs ",
    "currently installed in the cache (see '?list_germline_dbs')."
)

.stop_on_existing_cached_germline_db <- function(db_name)
{
    msg1 <- c("Germline db ", db_name, " is already installed ",
              "in igblastr's persistent cache.")
    msg3 <- c("Use 'overwrite=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(.LIST_GERMLINE_DB_TIP), "\n  ", wmsg(msg3))
}

.install_germline_db_succeeded <- function(db_name)
{
    message(wmsg("Germline db ", db_name, " successfully installed ",
                 "in igblastr's persistent cache.", margin=0L))
    message(wmsg("Call use_germline_db(\"", db_name, "\") to select ",
                 "it as the germline db to use with igblastn(). ",
                 .LIST_GERMLINE_DB_TIP, margin=0L))
}

.warn_if_auxdata_not_added_to_germline_db <- function(db_path)
{
    if (file.exists(make_germline_db_auxdata_path(db_path)))
        return()
    db_name <- basename(db_path)
    organism <- infer_organism_shortname_from_db_name(db_name)
    msg1 <- c("The auxiliary data didn't get added to ",
              "germline db ", db_name, " because the \"coding frame start\" ",
              "could not be determined for some of the J alleles in the db.")
    if (is.na(organism)) {
        msg2 <- "no auxiliary data"
    } else {
        msg2 <- c("the auxiliary data included in IgBLAST for ", organism)
    }
    msg2 <- c("This means that, when you select the db as the ",
              "germline db to use with igblastn(), ", msg2, " will ",
              "be used by default. ",
              "See documentation of the 'auxiliary_data' argument ",
              "in '?igblastn'.")
    warning(wmsg(msg1), "\n  ", wmsg(msg2))
}

### Not exported!
install_germline_db <- function(install_dir, db_name, fasta_dir, loci,
                                gapped=FALSE, with.intdata=FALSE,
                                excluded_J_alleles=NULL,
                                with.auxdata=FALSE, imgt.fasta=FALSE,
                                disambiguate.allele.names=FALSE,
                                if.exists=c("error", "overwrite", "no-op"),
                                verbose=FALSE, cheer.if.success=FALSE)
{
    stopifnot(isSingleNonWhiteString(install_dir), dir.exists(install_dir),
              isSingleNonWhiteString(db_name))
    if (!isSingleNonWhiteString(fasta_dir))
        stop(wmsg("'fasta_dir' must be a single (non-empty) string"))
    if (!dir.exists(fasta_dir))
        stop(wmsg(fasta_dir, ": directory not found"))
    checkarg_loci(loci)
    if (!isTRUEorFALSE(gapped))
        stop(wmsg("'gapped' must be TRUE or FALSE"))
    checkarg_with.intdata(with.intdata, gapped)
    if (!isTRUEorFALSE(with.auxdata))
        stop(wmsg("'with.auxdata' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(imgt.fasta))
        stop(wmsg("'imgt.fasta' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))
    if.exists <- match.arg(if.exists)
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
    stopifnot(isTRUEorFALSE(cheer.if.success))

    db_path <- file.path(install_dir, db_name)
    if (dir.exists(db_path)) {
        if (if.exists == "no-op")
            return()
        if (if.exists == "error")
            .stop_on_existing_cached_germline_db(db_name)
    }

    if (verbose) {
        what <- c("INSTALLING ", db_name)
        message("\n====== START ", wmsg(what), " ======\n")
    }

    create_germline_db(db_path, fasta_dir, loci,
                       gapped=gapped, with.intdata=with.intdata,
                       excluded_J_alleles=excluded_J_alleles,
                       with.auxdata=with.auxdata, imgt.fasta=imgt.fasta,
                       disambiguate.allele.names=disambiguate.allele.names,
                       overwrite=TRUE, verbose=verbose)

    if (verbose)
        message("====== DONE ", wmsg(what), " ======\n")

    if (cheer.if.success)
        .install_germline_db_succeeded(db_name)

    if (with.auxdata)
        .warn_if_auxdata_not_added_to_germline_db(db_path)

    invisible(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_custom_germline_db()
###

### A thin wrapper to install_germline_db().
install_custom_germline_db <- function(db_name, fasta_dir,
                                       tcr.db=FALSE, loci="auto",
                                       gapped=FALSE, with.intdata=FALSE,
                                       with.auxdata=FALSE, imgt.fasta=FALSE,
                                       disambiguate.allele.names=FALSE,
                                       overwrite=FALSE, verbose=FALSE)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    if (has_whitespace(db_name))
        stop(wmsg("'db_name' cannot contain whitespace characters"))
    if (!has_prefix(db_name, "cus"))
        stop(wmsg("'db_name' must start with \"cus\""))

    loci <- normalize_user_supplied_loci(loci, tcr.db=tcr.db,
                                         stop.if.missing.regions=TRUE)
    ## TODO: Issue warning if db_name does not contain loci prefix.

    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))

    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    if.exists <- if (overwrite) "overwrite" else "error"
    install_germline_db(germline_dbs_home, db_name, fasta_dir, loci,
                        gapped=gapped, with.intdata=with.intdata,
                        with.auxdata=with.auxdata, imgt.fasta=imgt.fasta,
                        disambiguate.allele.names=disambiguate.allele.names,
                        if.exists=if.exists, verbose=verbose)
}


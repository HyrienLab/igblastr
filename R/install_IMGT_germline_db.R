### =========================================================================
### install_IMGT_germline_db()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_IMGT_releases()
###

### Returns the IMGT/V-QUEST releases from newest to oldest (latest first).
list_IMGT_releases <- function(recache=FALSE)
{
    latest_release <- get_latest_IMGT_release(recache=recache)
    all_zips <- list_archived_IMGT_zips(recache=recache)
    archived_releases <- sub("^[^0-9]*([-0-9]+).*$", "\\1", all_zips)
    c(latest_release, sort(archived_releases, decreasing=TRUE))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_germline_db()
###

.normalize_IMGT_release <- function(release)
{
    if (!isSingleNonWhiteString(release))
        stop(wmsg("'release' must be a single (non-empty) string"))
    all_releases <- list_IMGT_releases()
    if (!(release %in% all_releases)) {
        stop(wmsg("\"", release, "\" is not a valid IMGT/V-QUEST release."),
             "\n  ",
             wmsg("Latest IMGT/V-QUEST release is \"", all_releases[[1L]],
                  "\" (recommended). Use list_IMGT_releases() to list ",
                  "all releases."))
    }
    release
}

.normalize_IMGT_organism <- function(organism)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    chartr(" ", "_", organism)
}

.get_IMGT_release_local_store <- function(release)
{
    igblastr_cache <- R_user_dir("igblastr", "cache")
    path <- file.path(igblastr_cache, "store", "IMGT-releases", release)
}

.form_IMGT_germline_db_name <- function(release, organism="Homo sapiens",
                                        db_type=c("IG", "TR", "IG-TR"))
{
    organism <- .normalize_IMGT_organism(organism)
    db_type <- match.arg(db_type)
    sprintf("IMGT-%s.%s.%s", release, organism, db_type)
}

### Requires Perl.
install_IMGT_germline_db <- function(release, organism="Homo sapiens",
                                     db_type=c("IG", "TR", "IG-TR"),
                                     force=FALSE, ...)
{
    ## Check arguments.
    if (missing(release)) {
        all_releases <- list_IMGT_releases()
        stop(wmsg("Argument 'release' is required and must be set ",
                  "to a valid IMGT/V-QUEST release."),
             "\n  ",
             wmsg("Latest IMGT/V-QUEST release is \"", all_releases[[1L]],
                  "\" (recommended). Use list_VQUEST_releases() to list ",
                  "all releases."))
    }

    release <- .normalize_IMGT_release(release)
    organism <- .normalize_IMGT_organism(organism)
    db_type <- match.arg(db_type)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))

    ## Download IMGT/V-QUEST release to local store if it's not there already.
    local_store <- .get_IMGT_release_local_store(release)
    if (!dir.exists(local_store))
        download_and_unzip_IMGT_release(release, local_store, ...)

    ## Compute 'organism_path' and 'db_name'.
    organism_path <- find_organism_in_IMGT_store(organism, local_store)
    organism <- basename(organism_path)
    db_name <- .form_IMGT_germline_db_name(release, organism, db_type)

    ## Create IMGT germline db.
    germline_dbs <- get_germline_dbs_path()
    db_path <- file.path(germline_dbs, db_name)
    make_IMGT_germline_db(organism_path, db_path, db_type, force=force)

    ## Success!
    message("Germline db ", db_name, " successfully installed.")
    message("Call use_germline_db(\"", db_name, "\") to select it")
    message("as the germline database to use with igblastn().")

    invisible(db_name)
}


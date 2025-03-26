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
### list_IMGT_organisms()
###

.stop_on_missing_release <- function()
{
    all_releases <- list_IMGT_releases()
    stop(wmsg("Argument 'release' is required and must be set ",
              "to a valid IMGT/V-QUEST release."),
         "\n  ",
         wmsg("Latest IMGT/V-QUEST release is \"", all_releases[[1L]],
              "\" (recommended). Use list_IMGT_releases() to list ",
              "all releases."))
}

.path_to_IMGT_local_store <- function(release=NULL)
{
    local_store <- igblastr_cache(IMGT_LOCAL_STORE)
    if (!is.null(release)) {
        stopifnot(isSingleNonWhiteString(release))
        local_store <- file.path(local_store, release)
    }
    local_store
}

.validate_IMGT_release <- function(release)
{
    if (!isSingleNonWhiteString(release))
        stop(wmsg("'release' must be a single (non-empty) string"))
    ## First we try offline validation by checking the IMGT local store.
    if (dir.exists(.path_to_IMGT_local_store(release)))
        return(release)
    ## Off-line validation above failed so we try online validation.
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

list_IMGT_organisms <- function(release)
{
    if (missing(release))
        .stop_on_missing_release()
    release <- .validate_IMGT_release(release)

    ## Download IMGT/V-QUEST release to local store if it's not there already.
    local_store <- .path_to_IMGT_local_store(release)
    if (!dir.exists(local_store))
        download_and_unzip_IMGT_release(release, local_store)
    list_organisms_in_IMGT_local_store(local_store)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_germline_db()
###

### Requires Perl.
install_IMGT_germline_db <- function(release, organism="Homo sapiens",
                                     force=FALSE, ...)
{
    ## Check arguments.
    if (missing(release))
        .stop_on_missing_release()
    release <- .validate_IMGT_release(release)
    organism <- normalize_IMGT_organism(organism)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))

    ## Download IMGT/V-QUEST release to local store if it's not there already.
    local_store <- .path_to_IMGT_local_store(release)
    if (!dir.exists(local_store))
        download_and_unzip_IMGT_release(release, local_store, ...)

    ## Compute 'organism_path' and 'db_name'.
    organism_path <- find_organism_in_IMGT_local_store(organism, local_store)
    organism <- basename(organism_path)
    db_name <- form_IMGT_germline_db_name(release, organism)

    ## Create IMGT germline db.
    IG_path <- file.path(organism_path, "IG")
    germline_dbs_path <- get_germline_dbs_path(TRUE)  # guaranteed to exist
    db_path <- file.path(germline_dbs_path, db_name)
    create_germline_db(IG_path, db_path, force=force)

    ## Success!
    message("Germline db ", db_name, " successfully installed.")
    message("Call use_germline_db(\"", db_name, "\") to select it")
    message("as the germline database to use with igblastn().")

    invisible(db_name)
}


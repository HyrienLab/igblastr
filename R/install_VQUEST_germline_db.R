### =========================================================================
### install_VQUEST_germline_db()
### -------------------------------------------------------------------------


install_VQUEST_release_in_local_store <-
    function(release="LATEST", force=FALSE, ...)
{
    local_store <- get_local_VQUEST_store(release)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (dir.exists(local_store) && !force) {
        release <- basename(local_store)
        message("VQUEST ", release, " is already in local store. ",
                "Use 'force=TRUE' to reinstall.")
        return(invisible(FALSE))
    }
    download_and_unzip_VQUEST_release(release, local_store, ...)
    invisible(TRUE)
}

.make_VQUEST_germline_db_name <-
    function(release="LATEST", organism="Homo sapiens",
             db_type=c("IG", "TR", "IG-TR"))
{
    release <- basename(get_local_VQUEST_store(release))
    organism <- normalize_VQUEST_organism(organism)
    db_type <- match.arg(db_type)
    sprintf("VQUEST-%s.%s.%s", release, organism, db_type)
}

install_VQUEST_germline_db <-
    function(release="LATEST", organism="Homo sapiens",
             db_type=c("IG", "TR", "IG-TR"), force=FALSE, ...)
{
    local_store <- get_local_VQUEST_store(release)
    organism <- normalize_VQUEST_organism(organism)
    db_type <- match.arg(db_type)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (!dir.exists(local_store))
        download_and_unzip_VQUEST_release(release, local_store, ...)
    organism_path <- find_organism_in_VQUEST_store(organism, local_store)
    organism <- basename(organism_path)
    db_name <- .make_VQUEST_germline_db_name(release, organism, db_type)
    germline_dbs <- get_germline_dbs_path()
    db_path <- file.path(germline_dbs, db_name)
    create_VQUEST_db(organism_path, db_path, db_type, force=force)
    invisible(db_path)
}


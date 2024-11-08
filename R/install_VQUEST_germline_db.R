### =========================================================================
### install_VQUEST_germline_db()
### -------------------------------------------------------------------------


.get_germline_dbs_path <- function()
{
    file.path(R_user_dir("igblastr", "cache"), "germline_dbs")
}


install_VQUEST_germline_db <-
    function(organism="Homo sapiens", subdir=c("IG", "TR", "both"),
             release="LATEST", force=FALSE, ...)
{
    organism <- normalize_VQUEST_organism(organism)
    subdir <- match.arg(subdir)
    release <- normalize_VQUEST_release(release)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))

    local_store <- get_local_VQUEST_store()
    local_store_release <- file.path(local_store, release)
    if (!dir.exists(local_store_release)) {
        exdir <- download_and_unzip_VQUEST_archived_zip(release, ...)
        stopifnot(identical(exdir, local_store_release))  # sanity check
    }
    organism_path <- find_organism_in_VQUEST_store(organism, release)
    organism_path
}


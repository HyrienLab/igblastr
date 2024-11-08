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
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    local_store <- get_local_VQUEST_store(release)
    if (!dir.exists(local_store))
        download_and_unzip_VQUEST_release(release, local_store, ...)
    organism_path <- find_organism_in_VQUEST_store(organism, local_store)
    organism_path
}


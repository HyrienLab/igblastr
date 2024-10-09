
### Will work offline if the required IMGT germline sequences are already
### in the cache.
make_IMGT_igblast_germline_db <-
    function(organism="Homo_sapiens", subdir=c("IG", "TR", "both"),
             destdir=".", quiet=FALSE)
{
    organism <- normalize_IMGT_VQUEST_organism(organism)
    subdir <- match.arg(subdir)
    orgdir_cache <- IMGT_VQUEST_orgdir_cache(organism)
    orgsubdir_cache <- file.path(orgdir_cache, subdir)
    if (!dir.exists(orgsubdir_cache))
        download_IMGT_germline_sequences_to_cache(organism, subdir=subdir,
                                                  quiet=quiet)

}


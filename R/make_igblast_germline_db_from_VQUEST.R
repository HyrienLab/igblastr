
.germline_dbs_root <- function()
    file.path(R_user_dir("igblastr", "data"), "germline_dbs")

.copy_all_VQUEST_sequences_to_one_file <- function(orgsubdir_cache,
                                                   segment_label, destfile)
{
    stop("work in progress")
}

### Execute instructions provided at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### to set up a germline database from IMGT VQUEST germline sequences.
### Will work offline if the required VQUEST germline sequences are
### already in the cache.
make_igblast_germline_db_from_VQUEST <-
    function(organism="Homo_sapiens", subdir=c("IG", "TR", "both"),
             destdir=".", quiet=FALSE)
{
    organism <- normalize_VQUEST_organism(organism)
    subdir <- match.arg(subdir)
    orgdir_cache <- VQUEST_orgdir_cache(organism)
    orgsubdir_cache <- file.path(orgdir_cache, subdir)
    if (!dir.exists(orgsubdir_cache))
        download_VQUEST_germline_sequences_to_cache(organism, subdir=subdir,
                                                    quiet=quiet)

    ## Combine all V, all D and all J sequences, respectively, into separate
    ## files.
    tempdir <- tempfile()
    dir.create(tempdir)
    V_file <- file.path(tempdir, "V_sequences.fasta")
    .copy_all_VQUEST_sequences_to_one_file(orgsubdir_cache, "V", V_file)
    D_file <- file.path(tempdir, "D_sequences.fasta")
    .copy_all_VQUEST_sequences_to_one_file(orgsubdir_cache, "D", D_file)
    J_file <- file.path(tempdir, "J_sequences.fasta")
    .copy_all_VQUEST_sequences_to_one_file(orgsubdir_cache, "J", J_file)

    ## Pass each of V_file, D_file, and J_file thru
    ## IGBLAST_ROOT/bin/edit_imgt_file.pl and IGBLAST_ROOT/bin/makeblastdb
    ## See https://ncbi.github.io/igblast/cook/How-to-set-up.html for the
    ## exact commands to use.
    ## The 30 resulting files should end up in 'dbdir'.
    dbname <- paste0("VQUEST_", organism, "_", subdir)
    dbdir <- file.path(.germline_dbs_root(), dbname)

}


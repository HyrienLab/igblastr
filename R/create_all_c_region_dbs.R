### =========================================================================
### create_all_c_region_dbs()
### -------------------------------------------------------------------------


get_c_region_dbs_path <- function()
{
    igblastr_cache <- R_user_dir("igblastr", "cache")
    file.path(igblastr_cache, "c_region_dbs")
}

### Called at load time to create all the C-region dbs.
create_all_c_region_dbs <- function(force=FALSE)
{
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    c_region_dbs <- get_c_region_dbs_path()
    if (!dir.exists(c_region_dbs))
        dir.create(c_region_dbs, recursive=TRUE)

    ## Install IMGT C-region_dbs.
    IMGT_c_region_dir <- system.file(package="igblastr", "extdata",
                                     "constant_regions", "IMGT",
                                     mustWork=TRUE)
    organism_paths <- list.dirs(IMGT_c_region_dir, recursive=FALSE)
    for (organism_path in organism_paths) {
        organism <- basename(organism_path)
        db_name <- paste0("IMGT.", organism)
        db_path <- file.path(c_region_dbs, db_name)
        create_IMGT_c_region_db(organism_path, db_path, force=force)
    }

    ## Any other C-region dbs to install?
}


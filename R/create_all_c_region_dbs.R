### =========================================================================
### create_all_c_region_dbs()
### -------------------------------------------------------------------------


### TODO: Maybe move this to list_c_region_dbs.R, or move
### compile_all_c_region_dbs() and clean_all_c_region_dbs() from
### list_c_region_dbs.R to here.

### Do NOT call this in .onLoad()! It relies on create_IMGT_c_region_db()
### which requires that Perl and a valid IgBLAST installation (for
### the 'edit_imgt_file.pl' script) are already available on the machine.
### However, none of these things are guaranteed to be available at load-time,
### especially if it's the first time that the package gets loaded on the user
### machine (e.g. right after installing the package from source).
create_all_c_region_dbs <- function(destdir, force=FALSE)
{
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))

    ## Create IMGT C-region dbs.
    IMGT_c_region_dir <- system.file(package="igblastr", "extdata",
                                     "constant_regions", "IMGT",
                                     mustWork=TRUE)
    organism_paths <- list.dirs(IMGT_c_region_dir, recursive=FALSE)
    for (organism_path in organism_paths) {
        organism <- basename(organism_path)
        db_name <- paste0("IMGT.", organism)
        db_path <- file.path(destdir, db_name)
        create_IMGT_c_region_db(organism_path, db_path, force=force)
    }

    ## Any other C-region dbs to create?
}


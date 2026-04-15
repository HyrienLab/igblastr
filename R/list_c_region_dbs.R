### =========================================================================
### list_c_region_dbs() and related
### -------------------------------------------------------------------------


.warn_about_out_of_sync_c_region_dbs <- function()
{
    msg <- c("Your C-region dbs are out-of-sync with your version of ",
             "igblastr. Please reset them with 'reset_c_region_dbs()'. ",
             "Note that this will remove any user-installed C-region db! ",
             "See '?reset_c_region_dbs' for more information.")
    warning(wmsg(msg))
}

### Starting with igblastr 1.0.12/1.1.12, create_region_db() and related
### drop repeated alleles ("repeated" here means alleles with identical
### names **and** identical sequences). See .drop_repeated_alleles() in
### R/create_region_db.R for more info.
### As a consequence, built-in C-region db _IMGT.mouse.IGH.202509 is now
### expected to contain 55 alleles instead of 56. We use this as an
### indication that the C-region dbs need to be reset.
.warn_if_c_region_dbs_need_reset <- function(c_region_dbs_home)
{
    stopifnot(isSingleNonWhiteString(c_region_dbs_home),
              dir.exists(c_region_dbs_home))
    db_name <- "_IMGT.mouse.IGH.202509"
    db_path <- file.path(c_region_dbs_home, db_name)
    ok <- dir.exists(db_path)
    if (ok) {
        db_fasta_file <- get_db_fasta_file(db_path, "C")
        ok <- file.exists(db_fasta_file) &&
              length(fasta.seqlengths(db_fasta_file)) == 55L
    }
    if (!ok)
        .warn_about_out_of_sync_c_region_dbs()
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_c_region_dbs_home()
###

### Not exported!
### Returns path to C_REGION_DBS cache compartment (see R/cache-utils.R for
### details about igblastr's cache organization).
### When 'init.path=TRUE':
### - if the path to return exists then no further action is performed;
### - if the path to return does NOT exist then it's created and populated
###   with the built-in C-region dbs.
### This means that the returned path is only guaranteed to exist
### when 'init.path' is set to TRUE.
get_c_region_dbs_home <- function(init.path=FALSE)
{
    stopifnot(isTRUEorFALSE(init.path))
    c_region_dbs_home <- igblastr_cache(C_REGION_DBS)
    if (init.path && !dir.exists(c_region_dbs_home))
        reset_c_region_dbs()
    c_region_dbs_home
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_c_region_dbs()
###

### 'long.listing' is ignored when 'names.only' is TRUE.
### Returns a c_region_dbs_df object (data.frame extension) by default.
list_c_region_dbs <- function(builtin.only=FALSE,
                              names.only=FALSE, long.listing=FALSE)
{
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    ans <- list_dbs(c_region_dbs_home, what="C-region",
                    builtin.only=builtin.only,
                    names.only=names.only, long.listing=long.listing)
    if (is.data.frame(ans))
        class(ans) <- c("c_region_dbs_df", class(ans))
    .warn_if_c_region_dbs_need_reset(c_region_dbs_home)
    ans
}

print.c_region_dbs_df <- function(x, ...)
{
    ## list_c_region_dbs() already called 'get_c_region_dbs_home(TRUE)' so
    ## it's reasonable to assume that 'c_region_dbs_home' already exists and
    ## was initialized. So we don't need to call 'get_c_region_dbs_home(TRUE)'
    ## again.
    c_region_dbs_home <- get_c_region_dbs_home()
    print_dbs_df(x, c_region_dbs_home, what="C-region")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_c_region_db_name()
###

.stop_on_invalid_c_region_db_name <- function(db_name)
{
    msg1 <- c("\"", db_name, "\" is not the name of a cached C-region db.")
    msg2 <- c("Use list_c_region_dbs() to list the C-region dbs ",
              "currently installed in the cache (see '?list_c_region_dbs').")
    stop(wmsg(msg1), "\n  ", wmsg(msg2))
}

### Not exported!
check_c_region_db_name <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    all_db_names <- list_c_region_dbs(names.only=TRUE)
    if (!(db_name %in% all_db_names))
        .stop_on_invalid_c_region_db_name(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_c_region_db_path()
###

### Not exported!
### Note that the returned path is NOT guaranteed to exist.
get_c_region_db_path <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    stopifnot(db_name != "USING")
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    file.path(c_region_dbs_home, db_name)             # NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rm_c_region_db()
###

rm_c_region_db <- function(db_name)
{
    check_c_region_db_name(db_name)
    if (has_prefix(db_name, "_"))
        stop(wmsg("cannot remove a built-in C-region db"))

    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    db_in_use_path <- get_db_in_use(c_region_dbs_home, what="C-region")
    if (db_in_use_path != "" && basename(db_in_use_path) == db_name)
        set_db_in_use("C-region", "")  # cancel current selection

    db_path <- get_c_region_db_path(db_name)
    nuke_file(db_path)
}


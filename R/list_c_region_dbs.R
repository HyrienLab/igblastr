### =========================================================================
### list_c_region_dbs() and related
### -------------------------------------------------------------------------


### Not exported!
### Creates and populates <igblastr-cached>/c_region_dbs/ if it doesn't
### exist yet. So the returned path is guaranteed to exist.
get_c_region_dbs_path <- function()
{
    igblastr_cache <- R_user_dir("igblastr", "cache")
    c_region_dbs <- file.path(igblastr_cache, "c_region_dbs")
    if (!dir.exists(c_region_dbs)) {
        dir.create(c_region_dbs, recursive=TRUE)
        create_all_c_region_dbs(c_region_dbs)
    }
    c_region_dbs
}

### Not exported!
get_c_region_db_path <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name), db_name != "USING")
    c_region_dbs <- get_c_region_dbs_path()  # path guaranteed to exist
    file.path(c_region_dbs, db_name)         # path NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_c_region_dbs()
###

### TODO: If 'names.only=FALSE', return a data.frame with one db per row that
### reports some details about each db (e.g. nb of regions in each group).
list_c_region_dbs <- function(names.only=FALSE)
{
    if (!isTRUEorFALSE(names.only))
        stop(wmsg("'names.only' must be TRUE or FALSE"))
    stopifnot(names.only)  # names.only=FALSE not ready yet!
    c_region_dbs <- get_c_region_dbs_path()  # path guaranteed to exist
    all_db_names <- setdiff(list.files(c_region_dbs), "USING")
    sort(all_db_names)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_c_region_db()
###

### Return the C-regions in a DNAStringSet object.
read_c_region_db <- function(db_name)
{
    stop("not reay yet")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compile_all_c_region_dbs()
### clean_all_c_region_dbs()
###
### TODO: Should we move this to create_all_c_region_dbs.R?

### Not exported!
compile_all_c_region_dbs <- function()
{
    all_db_names <- list_c_region_dbs(names.only=TRUE)
    for (db_name in all_db_names) {
        db_path <- get_c_region_db_path(db_name)
        compile_region_db(db_path)
    }
}

### Not exported!
clean_all_c_region_dbs <- function()
{
    all_db_names <- list_c_region_dbs(names.only=TRUE)
    for (db_name in all_db_names) {
        db_path <- get_c_region_db_path(db_name)
        clean_region_db(db_path)
    }
}


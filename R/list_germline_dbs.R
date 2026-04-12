### =========================================================================
### list_germline_dbs() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_germline_dbs_home()
###

### Not exported!
### Returns path to GERMLINE_DBS cache compartment (see R/cache-utils.R for
### details about igblastr's cache organization).
### When 'init.path=TRUE':
### - if the path to return exists then no further action is performed;
### - if the path to return does NOT exist then it's created and populated
###   with the built-in germline dbs.
### This means that the returned path is only guaranteed to exist
### when 'init.path' is set to TRUE.
get_germline_dbs_home <- function(init.path=FALSE)
{
    stopifnot(isTRUEorFALSE(init.path))
    germline_dbs_home <- igblastr_cache(GERMLINE_DBS)
    if (init.path) {
        if (dir.exists(germline_dbs_home)) {
            install_missing_builtin_germline_dbs(germline_dbs_home)
        } else {
            reset_germline_dbs()
        }
    }
    germline_dbs_home
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Handle stale built-in germline dbs
###

.warn_about_new_set_of_builtin_dbs <- function()
{
    msg1 <- c("The set of built-in germline dbs has changed in this new ",
              "version of igblastr.")
    msg2 <- "Here are the changes that happened:"
    msg3 <- c("1. Some human and mouse dbs have been removed because the ",
              "annotations provided by OGRDB for the V alleles in these ",
              "dbs are problematic.")
    msg4 <- c("2. The db for rhesus monkey was updated to use version 2 ",
              "of OGRDB germline sets IGH_VDJ, IGK_VJ, and IGL_VJ (released ",
              "on Feb 5, 2026).")
    msg5 <- c("3. All the db names are now prefixed with \"_OGRDB.\" instead ",
              "of \"_AIRR.\". This better reflects the provenance of the ",
              "data included in the dbs.")
    msg6 <- c("Note that the old _AIRR.* germline dbs were kept. You can ",
              "remove them individually with rm_germline_db(). This warning ",
              "will only go away when all the _AIRR.* dbs are removed.")
    msg7 <- c("Alternatively, you can call reset_germline_dbs() to remove ",
              "all the old _AIRR.* germline dbs at once. ",
              "WARNING: This will also remove any non built-in germline db ",
              "that you installed with install_IMGT_germline_db() or ",
              "install_custom_germline_db()! See '?reset_germline_dbs' for ",
              "more information.")
    old_warning_length <- getOption("warning.length")
    options(warning.length=1500)
    on.exit(options(warning.length=old_warning_length))
    warning(wmsg(msg1),
            "\n\n  ", wmsg(msg2),
            "\n    ", wmsg(msg3, margin=7),
            "\n    ", wmsg(msg4, margin=7),
            "\n    ", wmsg(msg5, margin=7),
            "\n\n  ", wmsg(msg6),
            "\n\n  ", wmsg(msg7))
}

### Not exported!
warn_if_old_AIRR_dbs_are_present <- function()
{
    germline_dbs_home <- get_germline_dbs_home()
    old_AIRR_dbs <- list.files(germline_dbs_home, pattern="^_AIRR\\.")
    if (length(old_AIRR_dbs) != 0L)
        .warn_about_new_set_of_builtin_dbs()
}

### Not exported!
### In igblastr 1.0.26 and 1.1.26, _OGRDB.mouse.NOD_ShiLtJ.IGH+IGK+IGL.202501
### was replaced with _OGRDB.mouse.NOD_ShiLtJ.IGH+IGK+IGL.202205 and
### _OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202501 was replaced with
### _OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410, so the original germline dbs
### can go away.
RENAMED_BUILTIN_GERMLINE_DBS <- c(
    `202205`="_OGRDB.mouse.NOD_ShiLtJ.IGH+IGK+IGL.202501",
    `202410`="_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202501"
)

### Not exported!
warn_about_renamed_builtin_germline_dbs <- function()
{
    msg1 <- c("Germline dbs ", RENAMED_BUILTIN_GERMLINE_DBS[[1L]], " and ",
              RENAMED_BUILTIN_GERMLINE_DBS[[2L]], " were renamed ",
              "_OGRDB.mouse.NOD_ShiLtJ.IGH+IGK+IGL.202205 and ",
              "_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410, respectively.")
    msg2 <- c("The new names better reflect the release date of the ",
              "OGRDB germline sets that they contain. Note that the two ",
              "original germline dbs were kept for backward compatibility.")
    msg3 <- c("To avoid this warning, use the new germline dbs instead ",
              "(they're identical to the old ones). You can use ",
              "rm_germline_db() to remove the original germline dbs.")
    warning(wmsg(msg1), "\n  ", wmsg(msg2), "\n\n  ", wmsg(msg3))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_germline_dbs()
###

### 'long.listing' is ignored when 'names.only' is TRUE.
### Returns a germline_dbs_df object (data.frame extension) by default.
list_germline_dbs <- function(builtin.only=FALSE,
                              with.intdata.only=FALSE,
                              with.auxdata.only=FALSE,
                              names.only=FALSE, long.listing=FALSE)
{
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    ans <- list_dbs(germline_dbs_home, what="germline",
                    builtin.only=builtin.only,
                    with.intdata.only=with.intdata.only,
                    with.auxdata.only=with.auxdata.only,
                    names.only=names.only, long.listing=long.listing)
    if (is.data.frame(ans))
        class(ans) <- c("germline_dbs_df", class(ans))
    if (!names.only)
        warn_if_old_AIRR_dbs_are_present()
    ans
}

print.germline_dbs_df <- function(x, ...)
{
    ## list_germline_dbs() already called 'get_germline_dbs_home(TRUE)' so
    ## it's reasonable to assume that 'germline_dbs_home' already exists and
    ## was initialized. So we don't need to call 'get_germline_dbs_home(TRUE)'
    ## again (it would call install_missing_builtin_germline_dbs() which has
    ## a small non-negligible cost).
    germline_dbs_home <- get_germline_dbs_home()
    print_dbs_df(x, germline_dbs_home, what="germline")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### valid_germline_db_name()
### check_germline_db_name()
###

.stop_on_no_installed_germline_db_yet <- function()
{
    msg <- c("You don't have any installed germline database yet. ",
             "Use any of the install_*_germline_db() function (e.g. ",
             "install_IMGT_germline_db()) to install at least one.")
    stop(wmsg(msg))
}

.list_germline_db_names <- function()
{
    all_db_names <- list_germline_dbs(names.only=TRUE)
    ## Should never happen because of the built-in OGRDB dbs but we keep this
    ## check anyways just in case we get rid of the built-in OGRDB dbs in the
    ## future.
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    all_db_names
}

### Not exported!
valid_germline_db_name <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name))
    db_name %in% .list_germline_db_names()
}

.stop_on_invalid_germline_db_name <- function(db_name)
{
    msg1 <- c("\"", db_name, "\" is not the name of a cached germline db.")
    msg2 <- c("Use list_germline_dbs() to list the germline dbs ",
              "currently installed in the cache (see '?list_germline_dbs').")
    msg3 <- c("Note that you can use any of the install_*_germline_db() ",
              "function (e.g. install_IMGT_germline_db()) to install ",
              "additional germline dbs in the cache.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### Not exported!
check_germline_db_name <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    if (!valid_germline_db_name(db_name))
        .stop_on_invalid_germline_db_name(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_germline_db_path()
###

### Not exported!
### Note that the returned path is NOT guaranteed to exist.
get_germline_db_path <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    stopifnot(db_name != "USING")
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    file.path(germline_dbs_home, db_name)             # NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rm_germline_db()
###

### Built-in germline dbs are not removable except for the stale ones.
.germline_db_is_removable <- function(db_name)
{
    if (!has_prefix(db_name, "_") || has_prefix(db_name, "_AIRR."))
        return(TRUE)
    db_name %in% RENAMED_BUILTIN_GERMLINE_DBS
}

rm_germline_db <- function(db_name)
{
    check_germline_db_name(db_name)
    if (!.germline_db_is_removable(db_name))
        stop(wmsg("cannot remove a built-in germline db"))

    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    db_in_use_path <- get_db_in_use(germline_dbs_home, what="germline")
    if (db_in_use_path != "" && basename(db_in_use_path) == db_name)
        set_db_in_use("germline", "")  # cancel current selection

    db_path <- get_germline_db_path(db_name)
    nuke_file(db_path)
}


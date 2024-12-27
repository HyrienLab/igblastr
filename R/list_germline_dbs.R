### =========================================================================
### list_germline_dbs() and related
### -------------------------------------------------------------------------


get_germline_dbs_path <- function()
{
    igblastr_cache <- R_user_dir("igblastr", "cache")
    file.path(igblastr_cache, "germline_dbs")
}

get_germline_db_path <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name), db_name != "USING")
    germline_dbs <- get_germline_dbs_path()
    file.path(germline_dbs, db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_germline_dbs()
###

.count_germline_regions <- function(all_db_names, group=c("V", "D", "J"))
{
    group <- match.arg(group)
    vapply(all_db_names,
        function(db_name) {
            db_path <- get_germline_db_path(db_name)
            fasta_file <- file.path(db_path, paste0(group, ".fasta"))
            length(fasta.seqlengths(fasta_file))
        }, integer(1), USE.NAMES=FALSE)
}

list_germline_dbs <- function(names.only=FALSE)
{
    if (!isTRUEorFALSE(names.only))
        stop(wmsg("'names.only' must be TRUE or FALSE"))
    germline_dbs <- get_germline_dbs_path()
    if (dir.exists(germline_dbs)) {
        all_db_names <- sort(setdiff(list.files(germline_dbs), "USING"))
    } else {
        all_db_names <- character(0)
    }
    if (names.only)
        return(all_db_names)
    nVregions <- .count_germline_regions(all_db_names, group="V")
    nDregions <- .count_germline_regions(all_db_names, group="D")
    nJregions <- .count_germline_regions(all_db_names, group="J")
    data.frame(db_name=all_db_names,
               nVregions=nVregions, nDregions=nDregions, nJregions=nJregions)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### use_germline_db()
###

.stop_on_no_installed_germline_db_yet <- function()
{
    msg <- c("You don't have any installed germline database yet. ",
             "Use any of the install_*_germline_db() function (e.g. ",
             "install_IMGT_germline_db()) to install at least one.")
    stop(wmsg(msg))
}

.stop_on_no_selected_germline_db_yet <- function()
{
    msg <- c("You haven't selected the germline database to use ",
             "with igblastn() yet. Please select one with ",
             "use_germline_db(\"<db_name>\"). ",
             "See '?use_germline_db' for more information.")
    stop(wmsg(msg))
}

.get_germline_db_in_use <- function()
{
    all_db_names <- list_germline_dbs(TRUE)
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    germline_dbs <- get_germline_dbs_path()
    db_path <- get_db_in_use(germline_dbs, what="germline")
    if (db_path == "")
        .stop_on_no_selected_germline_db_yet()
    make_blastdbs(db_path)
    basename(db_path)
}

.stop_on_invalid_germline_db_name <- function(db_name)
{
    msg1 <- c("\"", db_name, "\" is not the name of a cached germline db.")
    msg2 <- c("Use list_germline_dbs() to list the germline databases ",
              "currently installed in the cache (see '?list_germline_dbs').")
    msg3 <- c("Note that you can use any of the install_*_germline_db() ",
              "function (e.g. install_IMGT_germline_db()) to install ",
              "additional germline databases in the cache.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

use_germline_db <- function(db_name=NULL)
{
    if (is.null(db_name))
        return(.get_germline_db_in_use())

    ## Check 'db_name'.
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    all_db_names <- list_germline_dbs(TRUE)
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    if (!(db_name %in% all_db_names))
        .stop_on_invalid_germline_db_name(db_name)

    germline_dbs <- get_germline_dbs_path()
    db_path <- file.path(germline_dbs, db_name)
    make_blastdbs(db_path)

    using_path <- file.path(germline_dbs, "USING")
    writeLines(db_name, using_path)
    invisible(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_germline_blastdbs()
###

clean_germline_blastdbs <- function()
{
    all_db_names <- list_germline_dbs(TRUE)
    for (db_name in all_db_names) {
        db_path <- get_germline_db_path(db_name)
        clean_blastdbs(db_path)
    }
}


### =========================================================================
### list_germline_dbs() and related
### -------------------------------------------------------------------------


### Returns "<igblastr-cache>/germline_dbs".
### Note that the returned path is guaranteed to exist.
get_germline_dbs_path <- function()
{
    germline_dbs <- file.path(igblastr_cache(), "germline_dbs")
    if (!dir.exists(germline_dbs))
        dir.create(germline_dbs, recursive=TRUE)
    germline_dbs
}

### Note that the returned path is NOT guaranteed to exist.
get_germline_db_path <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name), db_name != "USING")
    germline_dbs <- get_germline_dbs_path()
    file.path(germline_dbs, db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_germline_dbs()
###

.count_germline_regions <- function(all_db_names, region_type=c("V", "D", "J"))
{
    region_type <- match.arg(region_type)
    vapply(all_db_names,
        function(db_name) {
            db_path <- get_germline_db_path(db_name)
            fasta_file <- file.path(db_path, paste0(region_type, ".fasta"))
            length(fasta.seqlengths(fasta_file))
        }, integer(1), USE.NAMES=FALSE)
}

list_germline_dbs <- function(names.only=FALSE)
{
    if (!isTRUEorFALSE(names.only))
        stop(wmsg("'names.only' must be TRUE or FALSE"))
    germline_dbs <- get_germline_dbs_path()
    all_db_names <- sort(setdiff(list.files(germline_dbs), "USING"))
    if (names.only)
        return(all_db_names)

    nVregions <- .count_germline_regions(all_db_names, region_type="V")
    nDregions <- .count_germline_regions(all_db_names, region_type="D")
    nJregions <- .count_germline_regions(all_db_names, region_type="J")

    used <- character(length(all_db_names))
    db_path <- get_db_in_use(germline_dbs, what="germline")
    if (db_path != "")
        used[all_db_names %in% basename(db_path)] <- "*"

    data.frame(db_name=all_db_names,
               nVregions=nVregions, nDregions=nDregions, nJregions=nJregions,
               used=used)
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
    msg2 <- c("Use list_germline_dbs() to list the germline dbs ",
              "currently installed in the cache (see '?list_germline_dbs').")
    msg3 <- c("Note that you can use any of the install_*_germline_db() ",
              "function (e.g. install_IMGT_germline_db()) to install ",
              "additional germline dbs in the cache.")
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
### load_germline_db()
###

.normarg_region_types <- function(region_types=NULL)
{
    if (is.null(region_types))
        return(c("V", "D", "J"))
    if (!is.character(region_types) || anyNA(region_types))
        stop(wmsg("'region_types' must be NULL or ",
                  "a character vector with no NAs"))
    region_types <- toupper(region_types)
    if (length(region_types) == 1L) {
        region_types <- safeExplode(region_types)
    } else if (any(nchar(region_types) != 1L)) {
        stop(wmsg("'region_types' must have single-letter elements"))
    }
    if (!all(region_types %in% c("V", "D", "J")))
        stop(wmsg("'region_types' can only contain letters V, D, or J"))
    region_types
}

### Returns the V, D, and/or J regions in a DNAStringSet object.
load_germline_db <- function(db_name, region_types=NULL)
{
    db_path <- get_germline_db_path(db_name)
    if (!dir.exists(db_path))
        .stop_on_invalid_germline_db_name(db_name)
    region_types <- .normarg_region_types(region_types)
    fasta_files <- vapply(region_types,
        function(type) file.path(db_path, paste0(type, ".fasta")),
        character(1), USE.NAMES=FALSE)
    readDNAStringSet(fasta_files)
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


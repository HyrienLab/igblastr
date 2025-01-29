### =========================================================================
### list_germline_dbs() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_builtin_germline_dbs()
###

### Do NOT call this in .onLoad()! It relies on create_germline_db()
### which requires that Perl and a valid IgBLAST installation (for
### the 'edit_imgt_file.pl' script) are already available on the machine.
### However, none of these things are guaranteed to be available at load-time,
### especially if it's the first time that the package gets loaded on the user
### machine (e.g. right after installing the package from source).
.create_builtin_germline_dbs <- function(destdir, force=FALSE)
{
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))

    ## Create AIRR germline db for Human.
    human_path <- system.file(package="igblastr", "extdata",
                              "germline_sequences", "AIRR", "human",
                              mustWork=TRUE)
    db_name <- form_AIRR_germline_db_name(human_path)
    db_path <- file.path(destdir, db_name)
    create_germline_db(human_path, db_path, force=force)

    ## Create AIRR germline dbs for Mouse strains.
    mouse_path <- system.file(package="igblastr", "extdata",
                              "germline_sequences", "AIRR", "mouse",
                              mustWork=TRUE)
    strain_paths <- list.dirs(mouse_path, recursive=FALSE)
    for (strain_path in strain_paths) {
        strain <- basename(strain_path)
        db_name <- form_AIRR_germline_db_name(mouse_path, strain=strain)
        db_path <- file.path(destdir, db_name)
        create_germline_db(strain_path, db_path, force=force)
    }

    ## Any other builtin germline dbs to create?
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_germline_dbs_path()
### reset_germline_dbs_cache()
### get_germline_db_path()
###

### Returns "<igblastr-cache>/germline_dbs".
### Note that the returned path is guaranteed to exist.
get_germline_dbs_path <- function(init.path=FALSE)
{
    germline_dbs <- file.path(igblastr_cache(), "germline_dbs")
    if (!dir.exists(germline_dbs) && init.path) {
        ## We first create the builtin germline dbs in a temporary folder,
        ## and, if successful, we replace 'germline_dbs' with the temporary
        ## folder. This achieves atomicity in case something goes wrong.
        tmp_germline_dbs <- tempfile("germline_dbs_")
        dir.create(tmp_germline_dbs, recursive=TRUE)
        on.exit(nuke_file(tmp_germline_dbs))
        .create_builtin_germline_dbs(tmp_germline_dbs)
        rename_file(tmp_germline_dbs, germline_dbs, replace=TRUE)
    }
    germline_dbs
}

### Returns the value returned by unlink().
reset_germline_dbs_cache <- function()
{
    germline_dbs <- get_germline_dbs_path()  # path NOT guaranteed to exist
    nuke_file(germline_dbs)
}

### Note that the returned path is NOT guaranteed to exist.
get_germline_db_path <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name), db_name != "USING")
    germline_dbs <- get_germline_dbs_path(TRUE)  # path guaranteed to exist
    file.path(germline_dbs, db_name)             # path NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_germline_dbs()
###

### Returns a named integer vector with GERMLINE_GROUPS as names.
.tabulate_germline_db_by_group <- function(db_name)
{
    db_path <- get_germline_db_path(db_name)
    region_types <- c("V", "D", "J")
    all_counts <- lapply(region_types,
        function(region_type) {
            fasta_file <- file.path(db_path, paste0(region_type, ".fasta"))
            seqids <- names(fasta.seqlengths(fasta_file))
            counts <- tabulate_germline_seqids_by_group(seqids)
            if (!all(has_suffix(names(counts)[counts != 0L], region_type)))
                warning(wmsg("some seq ids in '" , fasta_file, "' don't look ",
                             "like ", region_type, "-region germline seq ids"))
            counts
        })
    data <- unlist(all_counts, use.names=FALSE)
    m <- matrix(data, ncol=length(GERMLINE_GROUPS), byrow=TRUE,
                dimnames=list(region_types, GERMLINE_GROUPS))
    colSums(m)
}

### Returns a matrix with 1 row per germline db and 1 column per group.
.tabulate_germline_dbs_by_group <- function(db_names)
{
    all_counts <- lapply(db_names, .tabulate_germline_db_by_group)
    data <- unlist(all_counts, use.names=FALSE)
    if (is.null(data))
        data <- integer(0)
    matrix(data, ncol=length(GERMLINE_GROUPS), byrow=TRUE,
           dimnames=list(NULL, GERMLINE_GROUPS))
}

list_germline_dbs <- function(names.only=FALSE)
{
    if (!isTRUEorFALSE(names.only))
        stop(wmsg("'names.only' must be TRUE or FALSE"))
    germline_dbs <- get_germline_dbs_path(TRUE)  # path guaranteed to exist
    all_db_names <- sort(setdiff(list.files(germline_dbs), "USING"))
    if (names.only)
        return(all_db_names)
    basic_stats <- .tabulate_germline_dbs_by_group(all_db_names)
    ans <- data.frame(db_name=all_db_names, basic_stats)
    class(ans) <- c("germline_dbs_df", class(ans))
    ans
}

print.germline_dbs_df <- function(x, ...)
{
    germline_dbs <- get_germline_dbs_path(TRUE)  # path guaranteed to exist
    print_dbs_df(x, germline_dbs, what="germline")
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
    germline_dbs <- get_germline_dbs_path(TRUE)  # path guaranteed to exist
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

    germline_dbs <- get_germline_dbs_path()  # path guaranteed to exist
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

### Not exported!
clean_germline_blastdbs <- function()
{
    germline_dbs <- get_germline_dbs_path()  # path NOT guaranteed to exist
    if (dir.exists(germline_dbs)) {
        all_db_names <- list_germline_dbs(TRUE)
        for (db_name in all_db_names) {
            db_path <- get_germline_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}


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
.create_builtin_germline_dbs <- function(destdir)
{
    stopifnot(isSingleNonWhiteString(destdir), !dir.exists(destdir))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'destdir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_destdir <- tempfile("builtin_germline_dbs_")
    dir.create(tmp_destdir, recursive=TRUE)
    on.exit(nuke_file(tmp_destdir))

    ## Create AIRR germline db for Human.
    human_path <- system.file(package="igblastr",
                              "extdata", "germline_sequences", "AIRR", "human",
                              mustWork=TRUE)
    db_name <- form_AIRR_germline_db_name(human_path)
    db_path <- file.path(tmp_destdir, db_name)
    create_germline_db(human_path, db_path)

    ## Create AIRR germline dbs for Mouse strains.
    mouse_path <- system.file(package="igblastr",
                              "extdata", "germline_sequences", "AIRR", "mouse",
                              mustWork=TRUE)
    strain_paths <- list.dirs(mouse_path, recursive=FALSE)
    for (strain_path in strain_paths) {
        strain <- basename(strain_path)
        db_name <- form_AIRR_germline_db_name(mouse_path, strain=strain)
        db_path <- file.path(tmp_destdir, db_name)
        create_germline_db(strain_path, db_path)
    }

    ## Any other built-in germline dbs to create?

    ## Everyting went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_germline_dbs_path()
### get_germline_db_path()
###

### Returns path to GERMLINE_DBS cache compartment (see R/cache-utils.R for
### details about igblastr's cache organization).
### When 'init.path=TRUE':
### - if the path to return exists then no further action is performed;
### - if the path to return does NOT exist then it's created and populated
###   with the built-in germline dbs.
### This means that the returned path is only guaranteed to exist
### when 'init.path' is set to TRUE.
get_germline_dbs_path <- function(init.path=FALSE)
{
    stopifnot(isTRUEorFALSE(init.path))
    germline_dbs_path <- igblastr_cache(GERMLINE_DBS)
    if (!dir.exists(germline_dbs_path) && init.path)
        .create_builtin_germline_dbs(germline_dbs_path)
    germline_dbs_path
}

### Note that the returned path is NOT guaranteed to exist.
get_germline_db_path <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name), db_name != "USING")
    germline_dbs_path <- get_germline_dbs_path(TRUE)  # guaranteed to exist
    file.path(germline_dbs_path, db_name)             # NOT guaranteed to exist
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
            fasta_file <- get_db_fasta_file(db_path, region_type)
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

list_germline_dbs <- function(builtin.only=FALSE, names.only=FALSE)
{
    if (!isTRUEorFALSE(builtin.only))
        stop(wmsg("'builtin.only' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(names.only))
        stop(wmsg("'names.only' must be TRUE or FALSE"))
    germline_dbs_path <- get_germline_dbs_path(TRUE)  # guaranteed to exist
    all_db_names <- setdiff(list.files(germline_dbs_path), "USING")
    if (builtin.only)
        all_db_names <- all_db_names[has_prefix(all_db_names, "_")]
    all_db_names <- sort_db_names(all_db_names)
    if (names.only)
        return(all_db_names)
    basic_stats <- .tabulate_germline_dbs_by_group(all_db_names)
    ans <- data.frame(db_name=all_db_names, basic_stats)
    class(ans) <- c("germline_dbs_df", class(ans))
    ans
}

print.germline_dbs_df <- function(x, ...)
{
    germline_dbs_path <- get_germline_dbs_path(TRUE)  # guaranteed to exist
    print_dbs_df(x, germline_dbs_path, what="germline")
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
    all_db_names <- list_germline_dbs(names.only=TRUE)
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    germline_dbs_path <- get_germline_dbs_path(TRUE)  # guaranteed to exist
    db_path <- get_db_in_use(germline_dbs_path, what="germline")
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
    all_db_names <- list_germline_dbs(names.only=TRUE)
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    if (!(db_name %in% all_db_names))
        .stop_on_invalid_germline_db_name(db_name)

    germline_dbs_path <- get_germline_dbs_path()  # guaranteed to exist
    db_path <- file.path(germline_dbs_path, db_name)
    make_blastdbs(db_path)

    using_path <- file.path(germline_dbs_path, "USING")
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
        function(region_type) get_db_fasta_file(db_path, region_type),
        character(1), USE.NAMES=FALSE)
    readDNAStringSet(fasta_files)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_germline_blastdbs()
###

### Not exported!
clean_germline_blastdbs <- function()
{
    germline_dbs_path <- get_germline_dbs_path()  # NOT guaranteed to exist
    if (dir.exists(germline_dbs_path)) {
        all_db_names <- list_germline_dbs(names.only=TRUE)
        for (db_name in all_db_names) {
            db_path <- get_germline_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}


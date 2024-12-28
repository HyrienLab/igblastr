### =========================================================================
### list_c_region_dbs() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_builtin_c_region_dbs()
###

### Do NOT call this in .onLoad()! It relies on create_IMGT_c_region_db()
### which requires that Perl and a valid IgBLAST installation (for
### the 'edit_imgt_file.pl' script) are already available on the machine.
### However, none of these things are guaranteed to be available at load-time,
### especially if it's the first time that the package gets loaded on the user
### machine (e.g. right after installing the package from source).
.create_builtin_c_region_dbs <- function(destdir, force=FALSE)
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
        db_name <- paste0("IMGT.", organism, ".IG")
        db_path <- file.path(destdir, db_name)
        create_IMGT_c_region_db(organism_path, db_path, force=force)
    }

    ## Any other C-region dbs to create?
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_c_region_dbs_path()
### get_c_region_db_path()
###

### Returns "<igblastr-cache>/c_region_dbs".
### When 'init.path=TRUE', the returned path is created if it doesn't exist
### yet, and populated with the builtin C-region dbs.
### Note that the returned path is only guaranteed to exist when 'init.path'
### is set to TRUE.
.get_c_region_dbs_path <- function(init.path=FALSE)
{
    stopifnot(isTRUEorFALSE(init.path))
    igblastr_cache <- R_user_dir("igblastr", "cache")
    c_region_dbs <- file.path(igblastr_cache, "c_region_dbs")
    if (!dir.exists(c_region_dbs) && init.path) {
        dir.create(c_region_dbs, recursive=TRUE)
        .create_builtin_c_region_dbs(c_region_dbs)
    }
    c_region_dbs
}

### Note that the returned path is NOT guaranteed to exist.
### Not exported!
get_c_region_db_path <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name), db_name != "USING")
    c_region_dbs <- .get_c_region_dbs_path(TRUE)  # path guaranteed to exist
    file.path(c_region_dbs, db_name)              # path NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_c_region_dbs()
###

.count_c_regions <- function(all_db_names)
{
    vapply(all_db_names,
        function(db_name) {
            db_path <- get_c_region_db_path(db_name)
            fasta_file <- file.path(db_path, "C.fasta")
            length(fasta.seqlengths(fasta_file))
        }, integer(1), USE.NAMES=FALSE)
}

list_c_region_dbs <- function(names.only=FALSE)
{
    if (!isTRUEorFALSE(names.only))
        stop(wmsg("'names.only' must be TRUE or FALSE"))
    c_region_dbs <- .get_c_region_dbs_path(TRUE)  # path guaranteed to exist
    all_db_names <- sort(setdiff(list.files(c_region_dbs), "USING"))
    if (names.only)
        return(all_db_names)
    nCregions <- .count_c_regions(all_db_names)
    data.frame(db_name=all_db_names, nCregions=nCregions)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### use_c_region_db()
###

### Returns "" if no db is currently in use.
.get_c_region_db_in_use <- function()
{
    c_region_dbs <- .get_c_region_dbs_path(TRUE)  # path guaranteed to exist
    db_path <- get_db_in_use(c_region_dbs, what="C-region")
    if (db_path == "")
        return(db_path)
    make_blastdbs(db_path)
    basename(db_path)
}

.stop_on_invalid_c_region_db_name <- function(db_name)
{
    msg1 <- c("\"", db_name, "\" is not the name of a cached C-region db.")
    msg2 <- c("Use list_c_region_dbs() to list the C-region dbs ",
              "currently installed in the cache (see '?list_c_region_dbs').")
    stop(wmsg(msg1), "\n  ", wmsg(msg2))
}

### Passing 'db_name=""' will cancel the current selection.
use_c_region_db <- function(db_name=NULL)
{
    if (is.null(db_name))
        return(.get_c_region_db_in_use())

    ## Check 'db_name'.
    if (!isSingleString(db_name))
        stop(wmsg("'db_name' must be a single string"))

    if (db_name == "") {
        ## Cancel the current selection.
        c_region_dbs <- .get_c_region_dbs_path()  # path NOT guaranteed to exist
        if (dir.exists(c_region_dbs)) {
            using_path <- file.path(c_region_dbs, "USING")
            unlink(using_path)
        }
    } else {
        all_db_names <- list_c_region_dbs(names.only=TRUE)
        if (!(db_name %in% all_db_names))
            .stop_on_invalid_c_region_db_name(db_name)
        c_region_dbs <- .get_c_region_dbs_path()  # path guaranteed to exist
        db_path <- file.path(c_region_dbs, db_name)
        make_blastdbs(db_path)
        using_path <- file.path(c_region_dbs, "USING")
        writeLines(db_name, using_path)
    }
    invisible(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_c_region_db()
###

### Returns the C regions in a DNAStringSet object.
read_c_region_db <- function(db_name)
{
    db_path <- get_c_region_db_path(db_name)
    if (!dir.exists(db_path))
        .stop_on_invalid_c_region_db_name(db_name)
    fasta_file <- file.path(db_path, "C.fasta")
    readDNAStringSet(fasta_file)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_c_region_blastdbs()
###

### Not exported!
clean_c_region_blastdbs <- function()
{
    c_region_dbs <- .get_c_region_dbs_path()  # path NOT guaranteed to exist
    if (dir.exists(c_region_dbs)) {
        all_db_names <- list_c_region_dbs(names.only=TRUE)
        for (db_name in all_db_names) {
            db_path <- get_c_region_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}


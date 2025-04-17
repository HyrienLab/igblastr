### =========================================================================
### list_c_region_dbs() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_builtin_c_region_dbs()
###

### Do NOT call this in .onLoad()! It relies on create_c_region_db()
### which requires that Perl and a valid IgBLAST installation (for
### the 'edit_imgt_file.pl' script) are already available on the machine.
### However, none of these things are guaranteed to be available at load-time,
### especially if it's the first time that the package gets loaded on the user
### machine (e.g. right after installing the package from source).
.create_builtin_c_region_dbs <- function(destdir)
{
    stopifnot(isSingleNonWhiteString(destdir), !dir.exists(destdir))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'destdir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_destdir <- tempfile("builtin_c_region_dbs_")
    dir.create(tmp_destdir, recursive=TRUE)
    on.exit(nuke_file(tmp_destdir))

    ## Create IMGT C-region dbs.
    IMGT_c_region_dir <- system.file(package="igblastr",
                                     "extdata", "constant_regions", "IMGT",
                                     mustWork=TRUE)
    organism_paths <- list.dirs(IMGT_c_region_dir, recursive=FALSE)
    for (organism_path in organism_paths) {
        db_name <- form_IMGT_c_region_db_name(organism_path)
        db_path <- file.path(tmp_destdir, db_name)
        create_c_region_db(organism_path, db_path)
    }

    ## Any other built-in C-region dbs to create?

    ## Everyting went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_c_region_dbs_path()
### get_c_region_db_path()
###

### Returns path to C_REGION_DBS cache compartment (see R/cache-utils.R for
### details about igblastr's cache organization).
### When 'init.path=TRUE':
### - if the path to return exists then no further action is performed;
### - if the path to return does NOT exist then it's created and populated
###   with the built-in C-region dbs.
### This means that the returned path is only guaranteed to exist
### when 'init.path' is set to TRUE.
.get_c_region_dbs_path <- function(init.path=FALSE)
{
    stopifnot(isTRUEorFALSE(init.path))
    c_region_dbs_path <- igblastr_cache(C_REGION_DBS)
    if (!dir.exists(c_region_dbs_path) && init.path)
        .create_builtin_c_region_dbs(c_region_dbs_path)
    c_region_dbs_path
}

### Note that the returned path is NOT guaranteed to exist.
### Not exported!
get_c_region_db_path <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name), db_name != "USING")
    c_region_dbs_path <- .get_c_region_dbs_path(TRUE)  # guaranteed to exist
    file.path(c_region_dbs_path, db_name)              # NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_c_region_dbs()
###

### Returns a named integer vector with GENE_LOCI as names.
.tabulate_c_region_db_by_locus <- function(db_name)
{
    db_path <- get_c_region_db_path(db_name)
    fasta_file <- get_db_fasta_file(db_path, "C")
    seqids <- names(fasta.seqlengths(fasta_file))
    tabulate_c_region_seqids_by_locus(seqids)
}

### Returns a matrix with 1 row per C-region db and 1 column per locus.
.tabulate_c_region_dbs_by_locus <- function(db_names)
{
    all_counts <- lapply(db_names, .tabulate_c_region_db_by_locus)
    data <- unlist(all_counts, use.names=FALSE)
    if (is.null(data))
        data <- integer(0)
    matrix(data, ncol=length(GENE_LOCI), byrow=TRUE,
           dimnames=list(NULL, GENE_LOCI))
}

list_c_region_dbs <- function(builtin.only=FALSE, names.only=FALSE)
{
    if (!isTRUEorFALSE(builtin.only))
        stop(wmsg("'builtin.only' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(names.only))
        stop(wmsg("'names.only' must be TRUE or FALSE"))
    c_region_dbs_path <- .get_c_region_dbs_path(TRUE)  # guaranteed to exist
    all_db_names <- setdiff(list.files(c_region_dbs_path), "USING")
    if (builtin.only)
        all_db_names <- all_db_names[has_prefix(all_db_names, "_")]
    all_db_names <- sort_db_names(all_db_names)
    if (names.only)
        return(all_db_names)
    basic_stats <- .tabulate_c_region_dbs_by_locus(all_db_names)
    ans <- data.frame(db_name=all_db_names, basic_stats)
    class(ans) <- c("c_region_dbs_df", class(ans))
    ans
}

print.c_region_dbs_df <- function(x, ...)
{
    c_region_dbs_path <- .get_c_region_dbs_path(TRUE)  # guaranteed to exist
    print_dbs_df(x, c_region_dbs_path, what="C-region")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### use_c_region_db()
###

### Returns "" if no db is currently in use.
.get_c_region_db_in_use <- function()
{
    c_region_dbs_path <- .get_c_region_dbs_path(TRUE)  # guaranteed to exist
    db_path <- get_db_in_use(c_region_dbs_path, what="C-region")
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
        c_region_dbs_path <- .get_c_region_dbs_path()  # NOT guaranteed to exist
        if (dir.exists(c_region_dbs_path)) {
            using_path <- file.path(c_region_dbs_path, "USING")
            unlink(using_path)
        }
    } else {
        all_db_names <- list_c_region_dbs(names.only=TRUE)
        if (!(db_name %in% all_db_names))
            .stop_on_invalid_c_region_db_name(db_name)
        c_region_dbs_path <- .get_c_region_dbs_path()  # guaranteed to exist
        db_path <- file.path(c_region_dbs_path, db_name)
        make_blastdbs(db_path)
        using_path <- file.path(c_region_dbs_path, "USING")
        writeLines(db_name, using_path)
    }
    invisible(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_c_region_db()
###

### Returns the C regions in a DNAStringSet object.
load_c_region_db <- function(db_name)
{
    db_path <- get_c_region_db_path(db_name)
    if (!dir.exists(db_path))
        .stop_on_invalid_c_region_db_name(db_name)
    fasta_file <- get_db_fasta_file(db_path, "C")
    readDNAStringSet(fasta_file)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_c_region_blastdbs()
###

### Not exported!
clean_c_region_blastdbs <- function()
{
    c_region_dbs_path <- .get_c_region_dbs_path()  # NOT guaranteed to exist
    if (dir.exists(c_region_dbs_path)) {
        all_db_names <- list_c_region_dbs(names.only=TRUE)
        for (db_name in all_db_names) {
            db_path <- get_c_region_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}


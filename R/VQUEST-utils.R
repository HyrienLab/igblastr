### =========================================================================
### Low-level utilities to retrieve data from the VQUEST download site
### -------------------------------------------------------------------------
###
### The two main functions defined in this file are:
###   - download_VQUEST_germline_sequences()
###   - download_VQUEST_germline_sequences_to_cache()
###
### Both functions are exported. Nothing else is.
###


### Do not remove the trailing slash.
.VQUEST_DOWNLOAD_ROOT_URL <- "https://www.imgt.org/download/V-QUEST/"

### Do not remove the trailing slash.
.VQUEST_REFERENCE_DIRECTORY <-
    paste0(.VQUEST_DOWNLOAD_ROOT_URL, "IMGT_V-QUEST_reference_directory/")

.VQUEST_RELEASE_FILE <-
    paste0(.VQUEST_DOWNLOAD_ROOT_URL, "IMGT_vquest_release.txt")

### Do not remove the trailing slash.
.VQUEST_ARCHIVES_URL <-
    paste0(.VQUEST_DOWNLOAD_ROOT_URL, "archives/")

.IG_FILES <- paste0("IG",
    c("HV", "HD", "HJ", "KV", "KJ", "LV", "LJ"), ".fasta")

.TR_FILES <- paste0("TR",
    c("AV", "AJ", "BV", "BD", "BJ", "DV", "DD", "DJ", "GV", "GJ"), ".fasta")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_VQUEST_archived_zips()
###

.cached_VQUEST_archived_zips_listing <- new.env(parent=emptyenv())

.make_df_from_matrix_of_tds <- function(m)
{
    ## Drop empty columns.
    col_is_empty <- vapply(seq_len(ncol(m)),
                           function(j) all(is_white_str(m[ , j])),
                           logical(1))
    m <- m[ , !col_is_empty, drop=FALSE]

    ## Sanity checks.
    EXPECTED_COLNAMES <- c("Name", "Last modified", "Size")
    stopifnot(ncol(m) == 3L)
    stopifnot(identical(tolower(colnames(m)), tolower(EXPECTED_COLNAMES)))

    m <- m[has_suffix(m[ , 1L], ".zip"), , drop=FALSE]
    df <- as.data.frame(m)
    df[[2L]] <- as.Date(df[[2L]])
    df
}

### Returns a data.frame with 3 columns (Name, Last modified, Size)
### and 1 row per .zip file.
.list_VQUEST_archived_zips <- function()
{
    html <- getUrlContent(.VQUEST_ARCHIVES_URL, type="text", encoding="UTF-8")
    xml <- read_html(html)
    #listing <- html_text(html_elements(xml, "section table tr td a"))
    #listing[has_suffix(listing, ".zip")]
    all_ths <- html_text(html_elements(xml, "section table tr th"))
    all_tds <- html_text(html_elements(xml, "section table tr td"))
    EXPECTED_NCOL <- 5L
    m <- matrix(all_tds, ncol=EXPECTED_NCOL, byrow=TRUE)
    colnames(m) <- all_ths[seq_len(EXPECTED_NCOL)]
    .make_df_from_matrix_of_tds(m)
}

### If 'as.df' is TRUE then the listing is returned as a data.frame
### with 3 columns (Name, Last modified, Size) and 1 row per .zip file.
list_VQUEST_archived_zips <- function(as.df=FALSE, recache=FALSE)
{
    if (!isTRUEorFALSE(as.df))
        stop(wmsg("'as.df' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    listing <- .cached_VQUEST_archived_zips_listing[["LISTING"]]
    if (is.null(listing) || recache) {
        listing <- .list_VQUEST_archived_zips()
        .cached_VQUEST_archived_zips_listing[["LISTING"]] <- listing
    }
    if (!as.df)
        listing <- listing[ , 1L]
    listing
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions to support download_VQUEST_germline_sequences() and
### download_VQUEST_germline_sequences_to_cache()
###

### Unused at the moment.
.get_VQUEST_current_release <- function()
{
    content <- getUrlContent(.VQUEST_RELEASE_FILE)
    sub("^([^ ]*)(.*)$", "\\1", content)
}

normalize_VQUEST_organism <- function(organism)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    chartr(" ", "_", organism)
}

.get_VQUEST_orgdir_url <- function(organism)
{
    organism <- normalize_VQUEST_organism(organism)
    orgdir_url <- paste0(.VQUEST_REFERENCE_DIRECTORY, organism, "/")
    if (!urlExists(orgdir_url))
        stop(organism, ": no such organism\n  ",
             wmsg("See ", .VQUEST_REFERENCE_DIRECTORY, " for ",
                  "the list of supported organisms."))
    orgdir_url
}

.VQUEST_reference_directory_cache_root <- function()
    file.path(R_user_dir("igblastr", "cache"),
              "IMGT_V-QUEST_reference_directory")

VQUEST_orgdir_cache <- function(organism)
{
    organism <- normalize_VQUEST_organism(organism)
    file.path(.VQUEST_reference_directory_cache_root(), organism)
}

### 'orgdir_url', 'dry.run', and 'destdir' are trusted.
### In particular:
### - 'orgdir_url' is trusted to have the trailing slash, which
###   it should if it was obtained with .get_VQUEST_orgdir_url().
### - 'destdir' is trusted to exist.
.download_IG_or_TR_files <- function(orgdir_url, subdir=c("IG", "TR"),
                                     dry.run=FALSE, destdir=".",
                                     method, quiet=FALSE)
{
    subdir <- match.arg(subdir)
    files <- switch(subdir, "IG"=.IG_FILES, "TR"=.TR_FILES)
    subdir_url <- paste0(orgdir_url, subdir, "/")
    downloaded_files <- character(0)
    for (file in files) {
        file_url <- paste0(subdir_url, file)
        ## Not all organisms have all the IG or TR files.
        if (urlExists(file_url)) {
            if (!dry.run) {
                destfile <- file.path(destdir, file)
                download.file(file_url, destfile, method, quiet)
            }
            downloaded_files <- c(downloaded_files, file)
        }
    }
    downloaded_files
}

### Thin wrapper to .download_IG_or_TR_files() above that creates the IG
### or TR subdir inside 'destdir' and download the files to it.
### We're trying to achieve atomic behavior: in case of error during download,
### the IG or TR subdir is removed. So we get everything or nothing!
### However, note that in case of user interrupt (CTRL+C), we can end up
### with a partial download and corrupted files!
### FIXME: Make the behavior **really** atomic, under any circumstances.
### Like for .download_IG_or_TR_files() above, 'orgdir_url' and 'destdir'
### are trusted.
.download_IG_or_TR_subdir <- function(orgdir_url, subdir=c("IG", "TR"),
                                      destdir=".", method, quiet=FALSE)
{
    subdir <- match.arg(subdir)
    destsubdir <- file.path(destdir, subdir)
    unlink(destsubdir, recursive=TRUE)
    if (!suppressWarnings(dir.create(destsubdir)))
        stop(wmsg("failed to create '", destsubdir, "'"))
    files <- try(.download_IG_or_TR_files(orgdir_url, subdir=subdir,
                                          destdir=destsubdir,
                                          method=method, quiet=quiet))
    if (inherits(files, "try-error")) {
        unlink(destsubdir, recursive=TRUE)
        stop(as.character(files))
    }
    files
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_VQUEST_germline_sequences()
###

### Exported!
download_VQUEST_germline_sequences <-
    function(organism="Homo_sapiens", subdir=c("IG", "TR", "both"),
             dry.run=FALSE, destdir=".", method, quiet=FALSE)
{
    orgdir_url <- .get_VQUEST_orgdir_url(organism)
    subdir <- match.arg(subdir)
    if (!isTRUEorFALSE(dry.run))
        stop(wmsg("'dry.run' must be TRUE or FALSE"))
    if (!isSingleString(destdir))
        stop(wmsg("'destdir' must be a single string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (missing(method))
        method <- getOption("download.file.method", "auto")

    if (subdir == "both") {
        IG_files <- .download_IG_or_TR_files(orgdir_url, subdir="IG",
                                             dry.run=dry.run, destdir=destdir,
                                             method=method, quiet=quiet)
        TR_files <- .download_IG_or_TR_files(orgdir_url, subdir="TR",
                                             dry.run=dry.run, destdir=destdir,
                                             method=method, quiet=quiet)
        files <- c(IG_files, TR_files)
    } else {
        files <- .download_IG_or_TR_files(orgdir_url, subdir=subdir,
                                          dry.run=dry.run, destdir=destdir,
                                          method=method, quiet=quiet)
    }
    if (dry.run) files else invisible(files)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_VQUEST_germline_sequences_to_cache()
###

### Exported!
download_VQUEST_germline_sequences_to_cache <-
    function(organism="Homo_sapiens", subdir=c("IG", "TR", "both"),
             method, quiet=FALSE)
{
    organism <- normalize_VQUEST_organism(organism)
    orgdir_url <- .get_VQUEST_orgdir_url(organism)
    subdir <- match.arg(subdir)
    if (missing(method))
        method <- getOption("download.file.method", "auto")

    orgdir_cache <- VQUEST_orgdir_cache(organism)
    if (!dir.exists(orgdir_cache)) {
        dir.create(orgdir_cache, recursive=TRUE)
        destfile <- file.path(orgdir_cache, basename(.VQUEST_RELEASE_FILE))
        download.file(.VQUEST_RELEASE_FILE, destfile, method, quiet=quiet)
    }

    if (subdir == "both") {
        IG_files <- .download_IG_or_TR_subdir(orgdir_url, subdir="IG",
                                              destdir=orgdir_cache,
                                              method=method, quiet=quiet)
        TR_files <- .download_IG_or_TR_subdir(orgdir_url, subdir="TR",
                                              destdir=orgdir_cache,
                                              method=method, quiet=quiet)
        files <- c(IG_files, TR_files)
    } else {
        files <- .download_IG_or_TR_subdir(orgdir_url, subdir=subdir,
                                           destdir=orgdir_cache,
                                           method=method, quiet=quiet)
    }
    files
}


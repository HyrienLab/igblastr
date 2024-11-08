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
### .get_latest_VQUEST_release()
### .list_archived_VQUEST_zips()
### .list_archived_VQUEST_releases()
### .normalize_VQUEST_release()
###

.VQUEST_cache <- new.env(parent=emptyenv())

.fetch_latest_VQUEST_release <- function()
{
    content <- getUrlContent(.VQUEST_RELEASE_FILE)
    sub("^([^ ]*)(.*)$", "\\1", content)
}

.get_latest_VQUEST_release <- function(recache=FALSE)
{
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    release <- .VQUEST_cache[["LATEST_RELEASE"]]
    if (is.null(release) || recache) {
        release <- .fetch_latest_VQUEST_release()
        .VQUEST_cache[["LATEST_RELEASE"]] <- release
    }
    release
}

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
.fetch_list_of_archived_VQUEST_zips <- function()
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
.list_archived_VQUEST_zips <- function(as.df=FALSE, recache=FALSE)
{
    if (!isTRUEorFALSE(as.df))
        stop(wmsg("'as.df' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    listing <- .VQUEST_cache[["ARCHIVES_TABLE"]]
    if (is.null(listing) || recache) {
        listing <- .fetch_list_of_archived_VQUEST_zips()
        .VQUEST_cache[["ARCHIVES_TABLE"]] <- listing
    }
    if (!as.df)
        listing <- listing[ , 1L]
    listing
}

.list_archived_VQUEST_releases <- function(recache=FALSE)
{
    all_zips <- .list_archived_VQUEST_zips(recache=recache)
    sort(sub("^[^0-9]*([-0-9]+).*$", "\\1", all_zips), decreasing=TRUE)
}

.normalize_VQUEST_release <- function(release="LATEST")
{
    if (!isSingleNonWhiteString(release))
        stop(wmsg("'release' must be a single (non-empty) string"))
    if (release == "LATEST")
        return(release)
    archived_releases <- .list_archived_VQUEST_releases()
    if (release %in% archived_releases)
        return(release)
    all_in_1string <- paste0("\"", archived_releases, "\"", collapse=", ")
    stop(wmsg("'release' must be \"LATEST\" (recommended), or ",
              "one of the release numbers available at ",
              .VQUEST_ARCHIVES_URL, ", e.g. \"202405-2\"."),
         "\n  ",
         wmsg("All available releases: ", all_in_1string, "."),
         "\n  ",
         wmsg("Note that old releases have not been tested and ",
              "are not guaranteed to be compatible with the ",
              "igblastr package."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_local_VQUEST_store()
###

get_local_VQUEST_store <- function(release=NULL)
{
    path <- file.path(R_user_dir("igblastr", "cache"),
                      "store", "VQUEST-releases")
    if (!is.null(release)) {
        release <- .normalize_VQUEST_release(release)
        if (release == "LATEST")
            release <- .get_latest_VQUEST_release()
        path <- file.path(path, release)
    }
    path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_and_unzip_VQUEST_release()
###

.download_and_unzip_latest_VQUEST_zip <- function(exdir, ...)
{
    release <- .get_latest_VQUEST_release()
    refdir_zip_filename <- "IMGT_V-QUEST_reference_directory.zip"
    refdir_zip <- download_as_tempfile(.VQUEST_DOWNLOAD_ROOT_URL,
                                       refdir_zip_filename, ...)
    unlink(exdir, recursive=TRUE, force=TRUE)
    unzip(refdir_zip, exdir=exdir)
}

.get_archived_VQUEST_zip_from_release <- function(release)
{
    stopifnot(isSingleNonWhiteString(release))
    all_zips <- .list_archived_VQUEST_zips()
    idx <- grep(release, all_zips, fixed=TRUE)
    if (length(idx) == 0L)
        stop(wmsg("Anomaly: no .zip file found at ",
                  .VQUEST_ARCHIVES_URL, " for release ", release))
    if (length(idx) > 1L)
        stop(wmsg("Anomaly: more that one .zip file found at ",
                  .VQUEST_ARCHIVES_URL, " for release ", release))
    all_zips[[idx]]
}

.unzip_archived_VQUEST_zip <- function(zipfile, release, exdir)
{
    unlink(exdir, recursive=TRUE, force=TRUE)
    unzip(zipfile, exdir=exdir, junkpaths=TRUE)
    refdir_zip <- file.path(exdir, "IMGT_V-QUEST_reference_directory.zip")
    unzip(refdir_zip, exdir=exdir)
    unlink(refdir_zip)
}

.download_and_unzip_archived_VQUEST_zip <- function(release, exdir, ...)
{
    archived_zip_filename <- .get_archived_VQUEST_zip_from_release(release)
    archived_zipfile <- download_as_tempfile(.VQUEST_ARCHIVES_URL,
                                             archived_zip_filename, ...)
    .unzip_archived_VQUEST_zip(archived_zipfile, release, exdir)
}

### Download and unzip in 'exdir'.
download_and_unzip_VQUEST_release <- function(release, exdir, ...)
{
    if (dir.exists(exdir))
        unlink(exdir, recursive=TRUE, force=TRUE)
    if (release == "LATEST") {
        .download_and_unzip_latest_VQUEST_zip(exdir, ...)
    } else {
        .download_and_unzip_archived_VQUEST_zip(release, exdir, ...)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normalize_VQUEST_organism()
### find_organism_in_VQUEST_store()
###

normalize_VQUEST_organism <- function(organism)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    chartr(" ", "_", organism)
}

.list_organisms_in_VQUEST_store <- function(refdir)
{
    if (!dir.exists(refdir))
        stop(wmsg("Anomaly: directory ", refdir, " not found"))
    sort(list.files(refdir))
}

find_organism_in_VQUEST_store <- function(organism, local_store)
{
    refdir <- file.path(local_store, "IMGT_V-QUEST_reference_directory")
    all_organisms <- .list_organisms_in_VQUEST_store(refdir)
    idx <- match(tolower(organism), tolower(all_organisms))
    if (!is.na(idx))
        return(file.path(refdir, all_organisms[[idx]]))
    all_in_1string <- paste0("\"", all_organisms, "\"", collapse=", ")
    stop(wmsg(organism, ": organism not found in ",
              "VQUEST release ", basename(local_store), "."),
         "\n  ",
         wmsg("Available organisms: ", all_in_1string, "."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions to support download_VQUEST_germline_sequences() and
### download_VQUEST_germline_sequences_to_cache()
###

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


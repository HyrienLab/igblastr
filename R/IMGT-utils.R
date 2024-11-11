### =========================================================================
### Low-level utilities to retrieve data from the IMGT/V-QUEST download site
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


### Do not remove the trailing slash.
.VQUEST_DOWNLOAD_ROOT_URL <- "https://www.imgt.org/download/V-QUEST/"

### .VQUEST_REFERENCE_DIRECTORY
.VQUEST_REFERENCE_DIRECTORY <- "IMGT_V-QUEST_reference_directory"

.VQUEST_RELEASE_FILE_URL <-
    paste0(.VQUEST_DOWNLOAD_ROOT_URL, "IMGT_vquest_release.txt")

### Do not remove the trailing slash.
.VQUEST_ARCHIVES_URL <-
    paste0(.VQUEST_DOWNLOAD_ROOT_URL, "archives/")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_latest_IMGT_release()
### list_archived_IMGT_zips()
###

.IMGT_cache <- new.env(parent=emptyenv())

.fetch_latest_IMGT_release <- function()
{
    content <- getUrlContent(.VQUEST_RELEASE_FILE_URL)
    sub("^([^ ]*)(.*)$", "\\1", content)
}

get_latest_IMGT_release <- function(recache=FALSE)
{
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    release <- .IMGT_cache[["LATEST_RELEASE"]]
    if (is.null(release) || recache) {
        release <- .fetch_latest_IMGT_release()
        .IMGT_cache[["LATEST_RELEASE"]] <- release
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
.fetch_list_of_archived_IMGT_zips <- function()
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
list_archived_IMGT_zips <- function(as.df=FALSE, recache=FALSE)
{
    if (!isTRUEorFALSE(as.df))
        stop(wmsg("'as.df' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    listing <- .IMGT_cache[["ARCHIVES_TABLE"]]
    if (is.null(listing) || recache) {
        listing <- .fetch_list_of_archived_IMGT_zips()
        .IMGT_cache[["ARCHIVES_TABLE"]] <- listing
    }
    if (!as.df)
        listing <- listing[ , 1L]
    listing
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_and_unzip_IMGT_release()
###

.download_and_unzip_latest_IMGT_zip <- function(exdir, ...)
{
    release <- get_latest_IMGT_release()
    refdir_zip_filename <- paste0(.VQUEST_REFERENCE_DIRECTORY, ".zip")
    refdir_zip <- download_as_tempfile(.VQUEST_DOWNLOAD_ROOT_URL,
                                       refdir_zip_filename, ...)
    unlink(exdir, recursive=TRUE, force=TRUE)
    unzip(refdir_zip, exdir=exdir)
}

.get_archived_IMGT_zip <- function(release)
{
    stopifnot(isSingleNonWhiteString(release))
    all_zips <- list_archived_IMGT_zips()
    idx <- grep(release, all_zips, fixed=TRUE)
    if (length(idx) == 0L)
        stop(wmsg("Anomaly: no .zip file found at ",
                  .VQUEST_ARCHIVES_URL, " for release ", release))
    if (length(idx) > 1L)
        stop(wmsg("Anomaly: more that one .zip file found at ",
                  .VQUEST_ARCHIVES_URL, " for release ", release))
    all_zips[[idx]]
}

.unzip_archived_IMGT_zip <- function(zipfile, release, exdir)
{
    unlink(exdir, recursive=TRUE, force=TRUE)
    unzip(zipfile, exdir=exdir, junkpaths=TRUE)
    refdir_zip_filename <- paste0(.VQUEST_REFERENCE_DIRECTORY, ".zip")
    refdir_zip <- file.path(exdir, refdir_zip_filename)
    unzip(refdir_zip, exdir=exdir)
    unlink(refdir_zip)
}

.download_and_unzip_archived_IMGT_zip <- function(release, exdir, ...)
{
    archived_zip_filename <- .get_archived_IMGT_zip(release)
    archived_zipfile <- download_as_tempfile(.VQUEST_ARCHIVES_URL,
                                             archived_zip_filename, ...)
    .unzip_archived_IMGT_zip(archived_zipfile, release, exdir)
}

### Download and unzip in 'exdir'.
download_and_unzip_IMGT_release <- function(release, exdir, ...)
{
    if (dir.exists(exdir))
        unlink(exdir, recursive=TRUE, force=TRUE)
    if (release == get_latest_IMGT_release()) {
        .download_and_unzip_latest_IMGT_zip(exdir, ...)
    } else {
        .download_and_unzip_archived_IMGT_zip(release, exdir, ...)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_organism_in_IMGT_store()
###

.list_organisms_in_IMGT_store <- function(refdir)
{
    if (!dir.exists(refdir))
        stop(wmsg("Anomaly: directory ", refdir, " not found"))
    sort(list.files(refdir))
}

find_organism_in_IMGT_store <- function(organism, local_store)
{
    refdir <- file.path(local_store, .VQUEST_REFERENCE_DIRECTORY)
    all_organisms <- .list_organisms_in_IMGT_store(refdir)
    idx <- match(tolower(organism), tolower(all_organisms))
    if (!is.na(idx))
        return(file.path(refdir, all_organisms[[idx]]))
    all_in_1string <- paste0("\"", all_organisms, "\"", collapse=", ")
    stop(wmsg(organism, ": organism not found in ",
              "IMGT/V-QUEST release ", basename(local_store), "."),
         "\n  ",
         wmsg("Available organisms: ", all_in_1string, "."))
}


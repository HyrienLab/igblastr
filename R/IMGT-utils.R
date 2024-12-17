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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_IMGT_C_segments()
###

.VQUEST_C_SEGMENTS_URL <- "https://www.imgt.org/genedb/GENElect"

### All kinds of conventions are used across the IMGT website to name
### organisms. Picking one and sticking to it would boring I guess...
.map_organism_to_IMGT_species <- function(organism)
{
    stopifnot(isSingleNonWhiteString(organism))
    organism <- chartr("_", " ", organism)
    IMGT_species <- c("Homo sapiens", "Mus", "Rat",
                      "Vicugna pacos", "Oryctolagus cuniculus")
    m <- match(tolower(organism), tolower(IMGT_species))
    if (is.na(m)) {
        errmsg <- c("to the best of our knowledge, IMGT does not host ",
                    "C segments for organism ", organism)
        m <- switch(tolower(organism),
                    "human"=1L,
                    "mouse"=, "mus musculus"= 2L,
                    "rattus norvegicus"=3L,
                    "alpaca"=4L,
                    "rabbit"=5L,
                    stop(wmsg(errmsg))
        )
    }
    IMGT_species[m]
}

### 'fasta_lines' must be a character vector.
### Returns FALSE if no FASTA records or if all records are empty.
.is_dna_fasta <- function(fasta_lines)
{
    stopifnot(is.character(fasta_lines))
    header_idx <- grep("^>", fasta_lines)
    if (length(header_idx) == 0L)
        return(FALSE)  # no FASTA records
    dna_lines <- fasta_lines[- header_idx]
    dna <- paste(tolower(dna_lines), collapse="")
    if (!nzchar(dna))
        return(FALSE)  # all records are empty
    all(safeExplode(dna) %in% c("a", "c", "g", "t"))
}

### Fetch the C segments (as nucleotide sequences) from the links provided
### in the tables displayed at:
###   https://www.imgt.org/vquest/refseqh.html#constant-sets
### Unfortunately these links redirect us to ugly HTML pages that we need
### to scrape to extract the nucleotide sequences.
### Note that:
### - The IGHC group is available for Human, Mouse, Rat, Alpaca, and Rabbit.
### - The IGKC and IGLC groups are only available for Human.
### - The IMGT folks seem to use some kind of versioning convention for these
###   that is undocumented, unfortunately. As of Dec 16, 2024, these sequences
###   seem to be at version 14.1.
fetch_IMGT_C_segments <- function(organism, group=c("IGHC", "IGKC", "IGLC"),
                                  version="14.1")
{
    group <- match.arg(group)
    species <- .map_organism_to_IMGT_species(organism)
    stopifnot(isSingleNonWhiteString(version))
    query <- list(query=paste(version, group), species=species)
    html <- getUrlContent(.VQUEST_C_SEGMENTS_URL, query=query,
                          type="text", encoding="UTF-8")

    ## HTML document 'html' is expected to contain 2 <pre></pre> sections:
    ## - The first one is a section that describes the 15 fields of the
    ##   FASTA headers.
    ## - The second one contains our nucleotide sequences in FASTA format.
    ## To make the scrapping a little bit more robust we return the content
    ## of the first <pre></pre> section that contains valid FASTA.
    xml <- read_html(html)
    all_pres <- html_text(html_elements(xml, "pre"))
    for (pre in all_pres) {
        fasta_lines <- strsplit(pre, split="\n", fixed=TRUE)[[1L]]
        fasta_lines <- fasta_lines[nzchar(fasta_lines)]
        if (.is_dna_fasta(fasta_lines))
            return(fasta_lines)
    }
    constant_seq_url <- "https://www.imgt.org/vquest/refseqh.html#constant-sets"
    stop(wmsg("failed to fetch the nucleotide sequences of the C segments ",
              "for ", organism, " from the links provided in the tables ",
              "displayed at ", constant_seq_url))
}


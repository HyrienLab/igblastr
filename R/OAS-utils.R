### =========================================================================
### Some utilities to query OAS and manipulate OAS csv files
###
### - OAS csv files can be retrieved via the OAS web server at
###   https://opig.stats.ox.ac.uk/webapps/oas/
###
### - The OAS csv files are in MiAIRR format (Mi stands for “minimal”).
###   See https://onlinelibrary.wiley.com/doi/10.1002/pro.4205
### -------------------------------------------------------------------------
###


### This directory contains one subdirectory per study in the "paired"
### section of the OAS web site. (5 studies originally in Oct 2021, 12
### studies as of March 2025.)
.PAIRED_OAS_URL <- "https://opig.stats.ox.ac.uk/webapps/ngsdb/paired/"


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_OAS_csv_metadata()
### read_OAS_csv()
### extract_sequences_from_paired_OAS_df()
###

### Returns the metadata in a named list.
read_OAS_csv_metadata <- function(file)
{
    line <- readLines(file, n=1L)
    nc <- nchar(line)
    stopifnot(substr(line, 1L, 1L) == '"', substr(line, nc, nc) == '"')
    line <- substr(line, 2L, nc-1L)
    line <- gsub('""', '"', line, fixed=TRUE)
    jsonlite::fromJSON(line)
}

read_OAS_csv <- function(file, skip=1, ...)
{
    if (!isSingleNumber(skip))
        stop(wmsg("'skip' must be a single number"))
    tibble(read.table(file, header=TRUE, sep=",", skip=skip, ...))
}

.PAIRED_OAS_CORE_COLNAMES <- c("sequence_id_", "sequence_")
.suffixes <- rep(c("heavy", "light"), each=length(.PAIRED_OAS_CORE_COLNAMES))
.PAIRED_OAS_CORE_COLNAMES <- paste0(.PAIRED_OAS_CORE_COLNAMES, .suffixes)

.extract_paired_OAS_core_columns <- function(df)
{
    m <- match(.PAIRED_OAS_CORE_COLNAMES, colnames(df))
    missing_idx <- which(is.na(m))
    if (length(missing_idx) != 0L) {
        if ("sequence" %in% colnames(df))
            stop(wmsg("unpaired OAS format is not supported yet"))
        in1string <- paste0("\"", .PAIRED_OAS_CORE_COLNAMES[missing_idx], "\"",
                            collapse=", ")
        if (length(missing_idx) == 1L) {
            msg1 <- c("Column ", in1string, " is")
        } else {
            msg1 <- c("Columns ", in1string, " are")
        }
        msg1 <- c(msg1, " missing.")
        stop(wmsg(msg1, " Was this data.frame or tibble ",
                  "obtained with read_OAS_csv()?"))
    }
    df[ , m]
}

### We only support the "paired OAS" format for now, which is documented
### here: https://opig.stats.ox.ac.uk/webapps/oas/documentation_paired/
### Problem with the "unpaired OAS" format is that it doesn't include the
### sequence ids (required by igblastn()) so is not straightforward to support.
### Returns a DNAStringSet object of length 2 * nrow(df) where odd indices
### correspond to the heavy chain sequences and the even indices to the light
### chain sequences.
extract_sequences_from_paired_OAS_df <- function(df, add.prefix=FALSE)
{
    if (!is.data.frame(df))
        stop(wmsg("'df' must be a data.frame as obtained with read_OAS_csv()"))
    if (!isTRUEorFALSE(add.prefix))
        stop(wmsg("'add.prefix' must be TRUE or FALSE"))
    df <- .extract_paired_OAS_core_columns(df)
    ## Should never fail, unless OAS is doing something weird.
    stopifnot(anyDuplicated(df$sequence_id_heavy) == 0L,
              anyDuplicated(df$sequence_id_light) == 0L)
    sequences <- sequence_ids <- character(2L * nrow(df))
    sequences[c(TRUE, FALSE)] <- df$sequence_heavy
    sequences[c(FALSE, TRUE)] <- df$sequence_light
    sequence_ids[c(TRUE, FALSE)] <- df$sequence_id_heavy
    sequence_ids[c(FALSE, TRUE)] <- df$sequence_id_light
    ans <- DNAStringSet(setNames(sequences, sequence_ids))
    if (add.prefix) {
        ans_names <- paste0(c("heavy", "light"), "_chain_", names(ans))
        names(ans) <- ans_names
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_paired_OAS_studies()
### list_paired_OAS_units()
### download_paired_OAS_units()
###

.OAS_cache <- new.env(parent=emptyenv())

list_paired_OAS_studies <- function(as.df=FALSE, recache=FALSE)
{
    if (!isTRUEorFALSE(as.df))
        stop(wmsg("'as.df' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    listing <- .OAS_cache[["ALL_PAIRED_STUDIES"]]
    if (is.null(listing) || recache) {
        listing <- scrape_html_dir_index(.PAIRED_OAS_URL)
        .OAS_cache[["ALL_PAIRED_STUDIES"]] <- listing
    }
    if (!as.df)
        listing <- sub("/$", "", listing[ , "Name"])
    listing
}

### 'study' must be the name a subfolder in .PAIRED_OAS_URL e.g. "Jaffe_2022".
### If 'as.df' is TRUE then the listing is returned as a data.frame
### with 3 columns (Name, Last modified, Size) and 1 row per .csv.gz file.
list_paired_OAS_units <- function(study, as.df=FALSE, recache=FALSE)
{
    if (!isSingleNonWhiteString(study))
        stop(wmsg("'study' must be a single (non-empty) string"))
    if (!isTRUEorFALSE(as.df))
        stop(wmsg("'as.df' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    listing <- .OAS_cache[[study]]
    if (is.null(listing) || recache) {
        study_url <- paste0(.PAIRED_OAS_URL, study)
        csv_subdir <- if (study == "Phad_2022") "csv_paired" else "csv"
        url <- paste0(study_url, "/", csv_subdir, "/")
        listing <- scrape_html_dir_index(url, suffix=".csv.gz")
        .OAS_cache[[study]] <- listing
    }
    if (!as.df)
        listing <- listing[ , "Name"]
    listing
}

### Typical usage:
###   study <- "Jaffe_2022"
###   units <- list_paired_OAS_units(study)
###   download_paired_OAS_units(study, units[1:5])
download_paired_OAS_units <- function(study, units=NULL, destdir=".", ...)
{
    if (!isSingleNonWhiteString(study))
        stop(wmsg("'study' must be a single (non-empty) string"))
    if (is.null(units)) {
        units <- list_paired_OAS_units(study)
    } else {
        if (!is.character(units))
            stop(wmsg("'units' must be NULL or a character vector"))
    }
    for (unit in units) {
        url <- paste0(.PAIRED_OAS_URL, study, "/csv/", unit)
        destfile <- file.path(destdir, unit)
        download.file(url, destfile, ...)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_metadata_from_OAS_units()
### extract_sequences_from_paired_OAS_units()
###

### Returns the metadata of all CSV files in a data.frame with 1 row per file.
extract_metadata_from_OAS_units <- function(dir=".", pattern="\\.csv\\.gz$")
{
    files <- list.files(dir, pattern, full.names=TRUE)
    if (length(files) == 0L)
        stop(wmsg("no units to extract metadata from"))
    all_metadata <- lapply(files, read_OAS_csv_metadata)  # list of lists

    ## Trun 'all_metadata' into a list of 1-row data frames.
    all_metadata <- align_vectors_by_names(all_metadata)
    all_metadata <- lapply(all_metadata,
        function(metadata) {
            null_idx <- which(S4Vectors:::sapply_isNULL(metadata))
            metadata[null_idx] <- NA
            as.data.frame(metadata, optional=TRUE)
        })
    unit_names <- sub(pattern, "", basename(files))
    cbind(Unit=unit_names, do.call(rbind, all_metadata))
}

### Returns the sequences extracted from all CSV files in a DNAStringSet
### object. The names on the DNAStringSet object are the original sequence
### ids prefixed with the name of the unit that they are coming from.
extract_sequences_from_paired_OAS_units <-
    function(dir=".", pattern="\\.csv\\.gz$")
{
    files <- list.files(dir, pattern, full.names=TRUE)
    unit_names <- sub(pattern, "", basename(files))
    all_sequences <- lapply(setNames(files, unit_names),
        function(file) {
            df <- read_OAS_csv(file)
            extract_sequences_from_paired_OAS_df(df, add.prefix=TRUE)
        })
    ans <- do.call(c, unname(all_sequences))
    id_prefixes <- rep.int(unit_names, lengths(all_sequences))
    ans_ids <- paste0(id_prefixes, "_", names(ans))
    setNames(ans, ans_ids)
}


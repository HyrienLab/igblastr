### =========================================================================
### Low-level utilities to download data from OGRDB
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


OGRDB_URL <- "https://ogrdb.airr-community.org/"

### For some reason OGRDB URLs are expected to have the forward slash ("/")
### in their components encoded with "%25252f" instead of the usual "%2f".
### This means that we cannot simply use 'URLencode(x, reserved=TRUE)' to
### encode "/" in the components.
encode_OGRDB_URL_component <- function(component)
{
    stopifnot(isSingleString(component))
    gsub("/", "%25252f", URLencode(component), fixed=TRUE)
}

.OGRDB_organism_url <- function(organism, for.download=FALSE)
{
    stopifnot(isSingleNonWhiteString(organism), isTRUEorFALSE(for.download))
    top <- if (for.download) "download_germline_set" else "germline_sets"
    paste0(OGRDB_URL, top, "/", encode_OGRDB_URL_component(organism))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normalize_OGRDB_organism()
###

### Organisms available at OGRDB as of February 2026.
.OGRDB_ORGANISMS <- c(
    "Homo sapiens",
    "Macaca mulatta",
    "Mus musculus",
    "Oncorhynchus mykiss",
    "Peromyscus leucopus"
)

normalize_OGRDB_organism <- function(organism)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    organism <- chartr("_", " ", organism)
    m <- match(tolower(organism), tolower(.OGRDB_ORGANISMS))
    if (is.na(m)) {
        in1string <- paste0("\"", .OGRDB_ORGANISMS, "\"", collapse=", ")
        stop(wmsg("'organism' should be one of ", in1string))
    }
    .OGRDB_ORGANISMS[[m]]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Parse names of OGRDB germline sets
###

.IG_part_is_valid <- function(IG_part) grepl("^IG[HKL]", IG_part)

.stop_on_bad_OGRDB_set_name <- function(organism, set_name, what)
{
    org_url <- .OGRDB_organism_url(organism)
    msg1 <- c("Malformed OGRDB germline set name: can't extract ",
              "the ", what, " from ", organism, " germline set ",
              "name \"", set_name, "\".")
    msg2 <- c("Please visit ", org_url, " and make sure that the ",
              "name of the germline set is spelled correctly.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2))
}

### Returns a 2-column matrix with 1 row per set name. The column names
### are "strain" and "IG_part". The valid "IG_part" values for mouse are
### IGH, IGKV, IGKJ, IGLV, and IGLJ at the moment.
### Performs various sanity checks that should never fail for valid mouse
### germline set names.
.parse_OGRDB_mouse_set_names <- function(set_names)
{
    stopifnot(is.character(set_names), length(set_names) != 0L)

    stop_if_bad_set_names <- function(bad_idx) {
        if (length(bad_idx) != 0L) {
            idx1 <- bad_idx[[1L]]
            .stop_on_bad_OGRDB_set_name("Mus musculus", set_names[[idx1]],
                               "\"Species subgroup\" or \"IG* component\"")
        }
    }

    list_of_parts <- strsplit(set_names, " ", fixed=TRUE)
    nparts <- lengths(list_of_parts)
    bad_idx <- which(nparts < 2L | 3L < nparts)
    stop_if_bad_set_names(bad_idx)

    ## Special cases: "IGKJ (all strains)" and "IGLJ (all strains)".
    idx3 <- which(nparts == 3L)
    if (length(idx3) != 0L) {
        ok <- vapply(idx3,
                     function(i) identical(list_of_parts[[i]][2:3],
                                           c("(all", "strains)")),
                     logical(1))
        stop_if_bad_set_names(idx3[!ok])
        list_of_parts[idx3] <-
            lapply(idx3, function(i) c("", list_of_parts[[i]][[1L]]))
    }

    ## Construct 2-column matrix.
    ans <- matrix(unlist(list_of_parts), ncol=2L, byrow=TRUE,
                  dimnames=list(set_names, c("strain", "IG_part")))

    ## Check "IG_part" column.
    bad_idx <- which(!.IG_part_is_valid(ans[ , "IG_part"]))
    stop_if_bad_set_names(bad_idx)

    ans
}

### Note that OGRDB only defines the "Species subgroup" field for germline
### sets that are specific to a particular mouse strain at the moment, in
### which case the field is simply set to that strain. See:
###   https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
### When the "Species subgroup" field defined, it's always a **prefix** of
### the set name.
### Returns a character vector parallel to 'set_names' that is named with
### the supplied germline set names. The vector will contain empty strings
### for germline sets with an empty "Species subgroup" field.
extract_mouse_strains_from_OGRDB_set_names <- function(set_names)
{
    setNames(.parse_OGRDB_mouse_set_names(set_names)[ , "strain"], set_names)
}

### Returns a character vector parallel to 'set_names' that is named with
### the supplied germline set names.
### Note that **all** OGRDB germline set names have an "IG* component" that
### is guaranteed to start with a valid IG locus name, that is, with IGH,
### IGK, or IGL.
extract_IG_parts_from_OGRDB_set_names <- function(organism, set_names)
{
    stopifnot(isSingleNonWhiteString(organism), is.character(set_names))
    if (organism == "Mus musculus") {
        IG_parts <- .parse_OGRDB_mouse_set_names(set_names)[ , "IG_part"]
    } else {
        bad_idx <- which(!.IG_part_is_valid(set_names))
        if (length(bad_idx) != 0L) {
            set_name <- set_names[[bad_idx[[1L]]]]
            .stop_on_bad_OGRDB_set_name(organism, set_name, "\"IG*\" part")
        }
        IG_parts <- set_names
    }
    setNames(IG_parts, set_names)
}

### Returns a character vector parallel to 'set_names' that is named with
### the supplied germline set names.
extract_loci_from_OGRDB_set_names <- function(organism, set_names)
{
    substr(extract_IG_parts_from_OGRDB_set_names(organism, set_names), 1L, 3L)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .OGRDB_format2fileext()
### .normalize_OGRDB_format()
###

### Should always be called **before** .normalize_OGRDB_format() below.
.OGRDB_format2fileext <- function(format)
{
    stopifnot(isSingleNonWhiteString(format))
    switch(format, airr=".json", gapped=, ungapped=".fasta",
           stop(wmsg("unknown OGRDB format: ", format)))
}

.normalize_OGRDB_format <- function(format, organism, source_set=FALSE)
{
    if (!isSingleNonWhiteString(format))
        stop(wmsg("'format' must be a single (non-empty) string"))
    stopifnot(isSingleNonWhiteString(organism))
    if (!isTRUEorFALSE(source_set))
        stop(wmsg("'source_set' must be TRUE or FALSE"))

    if (source_set) {
        if (organism != "Homo sapiens")
            stop(wmsg("'source_set=TRUE' is only supported for Homo sapiens"))
    } else {
        if (organism == "Homo sapiens")
            format <- paste0(format, "_ex")
    }
    format
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_OGRDB_germline_set()
###

.build_OGRDB_germline_set_url <- function(organism, set_name, set_version,
                                          format)
{
    org_url <- .OGRDB_organism_url(organism, for.download=TRUE)
    stopifnot(isSingleNonWhiteString(set_name))
    if (organism == "Mus musculus") {
        strain <- extract_mouse_strains_from_OGRDB_set_names(set_name)
        if (strain != "")
            org_url <- paste0(org_url, "/", encode_OGRDB_URL_component(strain))
    }
    sprintf("%s/%s/%s/%s", org_url, encode_OGRDB_URL_component(set_name),
                           set_version, format)
}

### Note that the OGRDB website does NOT return an "HTTP 404 Not Found" error
### if the requested germline set does not exist. Instead it redirects to
### the website home page. However, looking at the headers of a HEAD requests
### gives us a clue.
.OGRDB_germline_set_exists <- function(url, ...)
{
    if (!has_internet())
        stop(wmsg("no internet"))
    config <- config(...)
    response <- try(HEAD(url, config, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(wmsg(as.character(response)))
    if (response$status_code == 404L)
        stop(wmsg("Not Found (HTTP 404): ", url))
    stop_for_status(response)
    headers <- headers(response)
    content_type <- headers[["content-type"]]
    content_disposition <- headers[["content-disposition"]]
    if (is.null(content_type) || is.null(content_disposition))
        return(FALSE)
    identical(content_type, "application/octet-stream")
}

### Similar to fetch_germline_set_from_OGRDB() (see R/OGRDB-API.R) but
### doesn't use OGRDB API.
### Returns the filename of the downloaded germline set.
download_OGRDB_germline_set <-
    function(organism, set_name, set_version,
             format=c("airr", "gapped", "ungapped"), source_set=FALSE,
             destdir=".", ...)
{
    ## Check arguments.
    organism <- normalize_OGRDB_organism(organism)
    if (!isSingleNonWhiteString(set_name))
        stop(wmsg("'set_name' must be a single (non-empty) string"))
    if (!isSingleNumber(set_version))
        stop(wmsg("'set_version' must be a single number"))
    if (!is.integer(set_version))
        set_version <- as.integer(set_version)
    format <- match.arg(format)
    fileext <- .OGRDB_format2fileext(format)
    format <- .normalize_OGRDB_format(format, organism, source_set=source_set)
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir)) {
        if (file.exists(destdir))
            stop(wmsg(destdir, ": not a directory"))
        stop(wmsg(destdir, ": no such directory"))
    }

    url <- .build_OGRDB_germline_set_url(organism, set_name, set_version,
                                         format)
    if (!.OGRDB_germline_set_exists(url)) {
        org_url <- .OGRDB_organism_url(organism)
        msg1 <- c(organism, " germline set \"", set_name, "\" ",
                  "(version ", set_version, ") not found.")
        msg2 <- c("Please visit ", org_url, " and make sure that the ",
                  "name of the germline set is spelled correctly and ",
                  "the specified version is valid.")
        stop(wmsg(msg1), "\n  ", wmsg(msg2))
    }
    filename <- paste0(format, fileext)
    download_file(url, file.path(destdir, filename), ...)
    filename
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_OGRDB_germline_set_to_OGRDB_store()
###

.encode_OGRDB_path_component <- function(path_component)
{
    stopifnot(isSingleNonWhiteString(path_component))
    gsub("/", ".slash.", chartr(" ", "_", path_component), fixed=TRUE)
}

.get_path_to_stored_OGRDB_germline_set <-
    function(organism, set_name, set_version, format, source_set=FALSE)
{
    organism <- normalize_OGRDB_organism(organism)
    if (!isSingleNonWhiteString(set_name))
        stop(wmsg("'set_name' must be a single (non-empty) string"))
    if (!isSingleNumber(set_version))
        stop(wmsg("'set_version' must be a single number"))
    if (!is.integer(set_version))
        set_version <- as.integer(set_version)
    fileext <- .OGRDB_format2fileext(format)
    format <- .normalize_OGRDB_format(format, organism, source_set=source_set)

    OGRDB_store <- igblastr_cache(OGRDB_STORE)
    path <- file.path(OGRDB_store, .encode_OGRDB_path_component(organism))
    if (organism == "Mus musculus") {
        strain <- extract_mouse_strains_from_OGRDB_set_names(set_name)
        if (strain != "")
            path <- file.path(path, .encode_OGRDB_path_component(strain))
    }
    set_name <- .encode_OGRDB_path_component(set_name)
    filename <- paste0(format, fileext)
    file.path(path, set_name, set_version, filename)
}

### Returns path to downloaded germline set.
download_OGRDB_germline_set_to_OGRDB_store <-
    function(organism, set_name, set_version,
             format=c("airr", "gapped", "ungapped"), source_set=FALSE,
             recache=FALSE, ...)
{
    format <- match.arg(format)
    local_file <- .get_path_to_stored_OGRDB_germline_set(organism,
                                            set_name, set_version,
                                            format, source_set=source_set)
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    ## Download OGRDB germline set to local store if it's not already there.
    if (!file.exists(local_file) || recache) {
        destdir <- dirname(local_file)
        if (!dir.exists(destdir)) {
            dir.create(destdir, recursive=TRUE)
            ## If the requested germline does not exist, the download below
            ## will fail and we will end up with an empty 'destdir'.
            on.exit(remove_empty_dir(destdir, parents=TRUE))
        }
        filename <- download_OGRDB_germline_set(organism,
                                   set_name, set_version,
                                   format=format, source_set=source_set,
                                   destdir=destdir, ...)
        ## Sanity check (should never fail).
        stopifnot(identical(filename, basename(local_file)))
    }
    local_file
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used by download_OGRDB_germline_sequences() (see
### R/download_OGRDB_germline_sequences.R) and
### download_OGRDB_germline_json() (see R/download_OGRDB_germline_json.R)
###

normalize_OGRDB_germline_sets <- function(germline_sets)
{
    if (!is.numeric(germline_sets) || length(germline_sets) == 0L)
        stop(wmsg("'germline_sets' must be a non-empty integer vector"))
    set_names <- names(germline_sets)
    if (is.null(set_names))
        stop(wmsg("'germline_sets' must have names"))
    if (anyNA(set_names))
        stop(wmsg("the names on 'germline_sets' cannot contain NAs"))
    set_names <- trimws2(set_names)
    if (any(nchar(set_names) == 0L))
        stop(wmsg("the names on 'germline_sets' cannot be empty"))
    if (anyDuplicated(set_names))
        stop(wmsg("the names on 'germline_sets' cannot contain duplicates"))
    if (!is.integer(germline_sets))
        germline_sets <- setNames(as.integer(germline_sets), set_names)
    if (anyNA(germline_sets))
        stop(wmsg("'germline_sets' cannot contain NAs"))
    germline_sets
}


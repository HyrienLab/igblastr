### =========================================================================
### Low-level utilities to retrieve data from the AIRR-community/OGRDB
### download site
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


### See short introduction to OGRDB REST API at
### https://wordpress.vdjbase.org/index.php/ogrdb_news/downloading-germline-sets-from-the-command-line-or-api/

.OGRDB_API_URL <- "https://ogrdb.airr-community.org/api"

.encode_OGRDB_API_query <- function(query)
{
    stopifnot(is.character(query))
    if (length(query) != 0L)
        stopifnot(!is.null(names(query)))
    vapply(query, function(x) gsub("/", "%25252f", URLencode(x)), character(1))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_germline_set_from_OGRDB()
###

### Does:
###   GET
###     /germline/set/
###       {species}/{species_subgroup}/{set_name}/{release_version}/{format}
### as documented at https://ogrdb.airr-community.org/api/
### Typical usage:
###   download_germline_set_from_OGRDB("Human", set_name="IGH_VDJ")
###   download_germline_set_from_OGRDB("Mouse", species_subgroup="C57BL/6",
###                                    set_name="C57BL/6 IGH")
download_germline_set_from_OGRDB <-
    function(species, species_subgroup=NULL, set_name,
             release_version="published",
             format=c("ungapped", "gapped", "airr",
                      "ungapped_ex", "gapped_ex", "airr_ex"))
{
    stopifnot(isSingleNonWhiteString(species),
              is.null(species_subgroup) ||
                  isSingleNonWhiteString(species_subgroup),
              isSingleNonWhiteString(set_name),
              isSingleNonWhiteString(release_version))
    format <- match.arg(format)

    query <- c(species=species, species_subgroup=species_subgroup,
               set_name=set_name, release_version=release_version,
               format=format)
    query <- .encode_OGRDB_API_query(query)
    url <- paste(c(.OGRDB_API_URL, "germline/set", query), collapse="/")
    content <- getUrlContent(url, type="text", encoding="UTF-8")
    ### Is it possible that 'content' will be NA? See R/REST_API.R in
    ### the UCSC.utils package for more info.
    stopifnot(isSingleString(content))

    if (!(format %in% c("airr", "airr_ex")))
        return(content)  # nucleotide sequences in FASTA format

    load_package_gracefully("jsonlite",
                            "setting 'format' to \"airr\" or \"airr_ex\"")
    parsed_json <- jsonlite::fromJSON(content, simplifyDataFrame=FALSE)
    ## Sanity checks.
    stopifnot(is.list(parsed_json),
              length(parsed_json) == 1L,
              identical(names(parsed_json), "GermlineSet"),
              is.list(parsed_json[[1L]]),
              length(parsed_json[[1L]]) == 1L)
    ans <- parsed_json[[1L]][[1L]]
    stopifnot(is.list(ans), !is.null(names(ans)),
              "allele_descriptions" %in% names(ans))
    ans$allele_descriptions <-
        jsonlite:::simplifyDataFrame(ans$allele_descriptions,
                                     flatten=FALSE, simplifyMatrix=TRUE)
    ans
}


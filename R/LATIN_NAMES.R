### A named character vector that represents a list of key/value pairs,
### where the keys are the short names used by IgBLAST for the 5 organisms
### that it officially supports. See function list_igblast_organisms() in
### file R/igblast_info.R.

LATIN_NAMES <- c(
    human="Homo sapiens",
    mouse="Mus musculus",
    rabbit="Oryctolagus cuniculus",
    rat="Rattus norvegicus",
    rhesus_monkey="Macaca mulatta"
)

### Treats 'organism' as a regular expression.
### Returns one of 'names(LATIN_NAMES)'.
lookup_igblast_organism <- function(organism)
{
    stopifnot(isSingleNonWhiteString(organism))
    ans <- grep(chartr(" ", "_", organism), names(LATIN_NAMES),
                ignore.case=TRUE, value=TRUE)
    if (length(ans) == 1L)
        return(ans)
    if (length(ans) >= 2L)
        stop(wmsg("ambiguous organism abbreviation: ", organism))
    idx <- grep(chartr("_", " ", organism), LATIN_NAMES,
                ignore.case=TRUE)
    if (length(idx) == 0L)
        stop(wmsg("unrecognized organism: ", organism))
    if (length(idx) >= 2L)
        stop(wmsg("ambiguous organism abbreviation: ", organism))
    names(LATIN_NAMES)[[idx]]
}

lookup_organism_latin_name <- function(organism)
    LATIN_NAMES[[lookup_igblast_organism(organism)]]

### Tries to map 'db_name' to one of the 5 organisms officially
### supported by IgBLAST.
infer_igblast_organism_from_db_name <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name))
    db_name <- chartr(".", "_", tolower(db_name))
    for (i in seq_along(LATIN_NAMES)) {
        igblast_organism <- names(LATIN_NAMES)[[i]]
        pattern <- paste0("_", chartr("_", ".", igblast_organism), "_")
        if (grepl(pattern, db_name, ignore.case=TRUE))
            return(igblast_organism)
        pattern <- paste0("_", chartr(" ", ".", LATIN_NAMES[[i]]), "_")
        if (grepl(pattern, db_name, ignore.case=TRUE))
            return(igblast_organism)
    }
    NA_character_
}


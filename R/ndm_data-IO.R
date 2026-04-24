### =========================================================================
### Read/write "ndm" files
### -------------------------------------------------------------------------
###


### Not exported!
### Not the true colnames used in IgBLAST "ndm" files: ours are all
### lowercase and we've replaced spaces with underscores.
### Note that many columns are redundant:
### - columns 'cdr1_start', 'fwr2_start', 'cdr2_start', and 'fwr3_start'
###   are redundant with columns 'fwr1_end', 'cdr1_end', 'fwr2_end',
###   and 'cdr2_end', respectively;
### - columns 'fwr1_start' and 'coding_frame_start' are redundant (and
###   column 'fwr1_start' is a dumb column anyways because it should always
###   be set to 1).
NDM_DATA_COL2CLASS <- c(
    allele_name="character",
    fwr1_start="integer",
    fwr1_end="integer",
    cdr1_start="integer",
    cdr1_end="integer",
    fwr2_start="integer",
    fwr2_end="integer",
    cdr2_start="integer",
    cdr2_end="integer",
    fwr3_start="integer",
    fwr3_end="integer",
    chain_type="character",
    coding_frame_start="integer"
)

### Not exported!
V_GENE_SEGMENTS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
V_GENE_DELINEATION_COLNAMES <- paste0(rep(V_GENE_SEGMENTS, each=2L),
                                      c("_start", "_end"))
stopifnot(all(V_GENE_DELINEATION_COLNAMES %in% names(NDM_DATA_COL2CLASS)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_ndm_data_col2class()
###

### Not exported!
check_ndm_data_col2class <- function(ndm_data, what="'ndm_data'")
{
    if (!is.data.frame(ndm_data))
        stop(wmsg(what, " must be a data.frame"))
    expected_colnames <- names(NDM_DATA_COL2CLASS)
    if (!identical(colnames(ndm_data), expected_colnames)) {
        in1string <- paste(expected_colnames, collapse=", ")
        stop(wmsg(what, " must have the following columns ",
                  "(in this order): ", in1string))
    }
    col2class <- vapply(ndm_data, function(x) class(x)[[1L]], character(1))
    if (!identical(col2class, NDM_DATA_COL2CLASS)) {
        in1string <- paste0("    ", names(NDM_DATA_COL2CLASS), " -> ",
                            NDM_DATA_COL2CLASS, collapse="\n")
        stop(wmsg(what, " must have the following column types:"), "\n",
             in1string)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_ndm_rows()
###

.validate_region_boundaries <- function(ndm_data, region, prev_end)
{
    starts <- ndm_data[ , paste0(region, "_start")]
    ends   <- ndm_data[ , paste0(region, "_end")]
    (starts == prev_end + 1L) & (ends > starts) & (ends %% 3L == 0L)
}

### Ignores the "chain_type" column at the moment.
validate_ndm_rows <- function(ndm_data, allow.repeated.rows=FALSE)
{
    check_ndm_data_col2class(ndm_data)
    if (!isTRUEorFALSE(allow.repeated.rows))
        stop(wmsg("'allow.repeated.rows' must be TRUE or FALSE"))
    if (allow.repeated.rows) {
        ## We allow duplicated entries in 'ndm_data' as long as they
        ## tell the same story.
        if (!rows_with_same_key_are_identical(ndm_data, "allele_name"))
            stop(wmsg("rows in 'ndm_data' with same \"allele_name\" ",
                      "must be identical"))
    } else {
        if (anyDuplicated(ndm_data[ , "allele_name"]))
            stop(wmsg("'ndm_data$allele_name' cannot contain duplicates"))
    }

    fwr1_ok <- .validate_region_boundaries(ndm_data, "fwr1", 0L)
    cdr1_ok <- .validate_region_boundaries(ndm_data, "cdr1", ndm_data$fwr1_end)
    fwr2_ok <- .validate_region_boundaries(ndm_data, "fwr2", ndm_data$cdr1_end)
    cdr2_ok <- .validate_region_boundaries(ndm_data, "cdr2", ndm_data$fwr2_end)
    fwr3_ok <- .validate_region_boundaries(ndm_data, "fwr3", ndm_data$cdr2_end)
    coding_frame_start_ok <- ndm_data[ , "coding_frame_start"] == 0L
    fwr1_ok & cdr1_ok & fwr2_ok & cdr2_ok & fwr3_ok & coding_frame_start_ok
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_ndm_data()
### write_ndm_data()
###

read_ndm_data <- function(filepath)
{
    read_broken_table(filepath, NDM_DATA_COL2CLASS)
}

write_ndm_data <- function(ndm_data, file="", check.data=FALSE)
{
    if (!isTRUEorFALSE(check.data))
        stop(wmsg("'check.data' must be TRUE or FALSE"))
    check_ndm_data_col2class(ndm_data)
    if (check.data) {
        ok <- validate_ndm_rows(ndm_data)
        if (!all(ok))
            stop(wmsg("'ndm_data' contains invalid rows. ",
                      "Use 'validate_ndm_rows()' to identify them."))
    }
    header <- paste0("#", paste(colnames(ndm_data), collapse=", "))
    cat(header, "\n", sep="", file=file)
    write.table(ndm_data, file, append=TRUE, quote=FALSE,
                sep="\t", row.names=FALSE, col.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### same_alleles_annot()
###

### Each of 'annot1' and 'annot2' must be a data.frame with an "allele_name"
### column. Typically used to compare two "intdata" or two "auxdata"
### data.frames.
### Returns TRUE if the 2 data.frames contain the same data (possibly
### with their rows in different order).
same_alleles_annot <- function(annot1, annot2)
{
    have_same_rows(annot1, annot2, "allele_name")
}


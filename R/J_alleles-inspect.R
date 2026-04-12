### =========================================================================
### Basic inspection of J allele sequences
### -------------------------------------------------------------------------
###


.get_auxdata_col <- function(auxdata, colname)
{
    if (!is.data.frame(auxdata))
        stop(wmsg("'auxdata' must be a data.frame as returned ",
                  "by load_auxdata() or compute_auxdata()"))
    if (!isSingleNonWhiteString(colname))
        stop(wmsg("'colname' must be a single (non-empty) string"))
    auxdata_col <- auxdata[[colname]]
    if (is.null(auxdata_col))
        stop(wmsg("'auxdata' has no \"", colname, "\" column. Make sure ",
                  "that it's a data.frame as returned by compute_auxdata(), ",
                  "compute_germline_db_auxdata(), or load_auxdata()."))
    auxdata_col
}

### Extracts the specified column from the 'auxdata' data.frame, and
### subset/reorder it to keep only the column values that correspond
### to the alleles in 'J_alleles'. Returns them in a named vector that
### is parallel to 'J_alleles' and has the allele names on it.
### The returned vector will have NAs for alleles that are not annotated
### in 'auxdata' or when 'auxdata[[colname]]' reports an NA for the allele.
.query_auxdata <- function(auxdata, J_alleles, colname)
{
    allele_names <- .get_auxdata_col(auxdata, "allele_name")
    if (!is(J_alleles, "DNAStringSet"))
        stop(wmsg("'J_alleles' must be DNAStringSet object"))
    J_names <- names(J_alleles)
    if (is.null(J_names))
        stop(wmsg("'J_alleles' must have names"))
    auxdata_col <- .get_auxdata_col(auxdata, colname)
    setNames(auxdata_col[match(J_names, allele_names)], J_names)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_J_alleles()
###

### Translates the coding frame contained in the J allele sequence.
### Only needs access to the "coding_frame_start" column in 'auxdata'.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'J_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'auxdata' or for which 'auxdata$coding_frame_start' has an NA.
translate_J_alleles <- function(J_alleles, auxdata)
{
    coding_frame_start <- .query_auxdata(auxdata, J_alleles,
                                         "coding_frame_start")
    ans <- rep.int(NA_character_, length(J_alleles))
    selection_idx <- which(!is.na(coding_frame_start))
    if (length(selection_idx) != 0L) {
        dna <- J_alleles[selection_idx]
        offset <- coding_frame_start[selection_idx]
        aa <- translate_codons(dna, offset=offset)
        ans[selection_idx] <- as.character(aa)
    }
    setNames(ans, names(J_alleles))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### J_allele_has_stop_codon()
###

### Only needs access to the "coding_frame_start" column in 'auxdata'.
### Returns a named logical vector that is parallel to 'J_alleles' and has
### the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'auxdata' or for which 'auxdata$coding_frame_start' has an NA.
J_allele_has_stop_codon <- function(J_alleles, auxdata)
{
    J_aa <- translate_J_alleles(J_alleles, auxdata)
    ans <- setNames(grepl("*", J_aa, fixed=TRUE), names(J_aa))
    ans[is.na(J_aa)] <- NA
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_fwr4()
###

### Only needs access to the "cdr3_end" column of the 'auxdata' data.frame.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'J_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'auxdata' or for which 'auxdata$cdr3_end' has an NA.
translate_fwr4 <- function(J_alleles, auxdata, max.codons=NA)
{
    if (!isSingleNumberOrNA(max.codons))
        stop(wmsg("'max.codons' must be a single number or NA"))
    if (!is.integer(max.codons))
        max.codons <- as.integer(max.codons)

    cdr3_end <- .query_auxdata(auxdata, J_alleles, "cdr3_end")  # 0-based
    ans <- rep.int(NA_character_, length(J_alleles))
    selection_idx <- which(!is.na(cdr3_end))
    if (length(selection_idx) != 0L) {
        dna <- J_alleles[selection_idx]
        offset <- cdr3_end[selection_idx] + 1L  # 0-based FWR4 start
        aa <- translate_codons(dna, offset=offset)
        ans[selection_idx] <- as.character(aa)
    }
    if (!is.na(max.codons))
        ans <- substr(ans, 1L, max.codons)
    setNames(ans, names(J_alleles))
}


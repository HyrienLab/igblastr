### =========================================================================
### Basic inspection of V allele sequences
### -------------------------------------------------------------------------
###


.INTDATA_AS_RETURNED_BY <- "load_intdata()"

### Not exported!
get_intdata_col <- function(intdata, colname)
{
    get_igdata_col(intdata, colname,
                   what="intdata", as_returned_by=.INTDATA_AS_RETURNED_BY)
}

### See query_igdata() in R/igdata-utils.R for what gets returned.
query_intdata <- function(intdata, V_alleles, colname, no.NAs=FALSE)
{
    if (!is(V_alleles, "DNAStringSet"))
        stop(wmsg("'V_alleles' must be a DNAStringSet object"))
    V_names <- names(V_alleles)
    if (is.null(V_names))
        stop(wmsg("'V_alleles' must have names"))
    query_igdata(intdata, V_names, colname, no.NAs=no.NAs,
                 what="intdata", as_returned_by=.INTDATA_AS_RETURNED_BY)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_V_alleles()
###

.translate_V_codons <- function(V_alleles, offsets, with.init.codon)
{
    stopifnot(is(V_alleles, "DNAStringSet"), is.integer(offsets),
              length(V_alleles) == length(offsets))
    ans <- rep.int(NA_character_, length(V_alleles))
    selection_idx <- which(!is.na(offsets))
    if (length(selection_idx) != 0L) {
        dna <- V_alleles[selection_idx]
        off <- offsets[selection_idx]
        aa <- translate_codons(dna, offset=off, with.init.codon=with.init.codon)
        ans[selection_idx] <- as.character(aa)
    }
    setNames(ans, names(V_alleles))
}

### Translates the coding frame contained in the V allele sequence.
### Only needs access to the "coding_frame_start" column in 'intdata'.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'V_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'intdata' or for which 'intdata$coding_frame_start' has
### an NA.
.translate_V_coding_frame <- function(V_alleles, intdata)
{
    offsets <- query_intdata(intdata, V_alleles, "coding_frame_start")
    .translate_V_codons(V_alleles, offsets, with.init.codon=TRUE)
}

### Only needs access to the "<V_segment>_start" and "<V_segment>_end"
### columns of the 'intdata' data.frame.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'V_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is
### not annotated in 'intdata' or for which 'intdata$<V_segment>_start'
### or 'intdata$<V_segment>_end' has an NA.
.translate_V_segment <- function(V_alleles, intdata, V_segment)
{
    .check_V_segment(V_segment)
    start_colname <- paste0(V_segment, "_start")
    end_colname <- paste0(V_segment, "_end")
    starts <- query_intdata(intdata, V_alleles, start_colname)  # 1-based
    ends <- query_intdata(intdata, V_alleles, end_colname)  # 1-based
    offsets <- starts - 1L
    with.init.codon <- V_segment == "fwr1"
    ans <- .translate_V_codons(V_alleles, offsets, with.init.codon)
    ncodons <- (ends - offsets) %/% 3L
    substr(ans, 1L, ncodons)
}

translate_V_alleles <- function(V_alleles, intdata, V_segment=NULL)
{
    if (is.null(V_segment))
        return(.translate_V_coding_frame(V_alleles, intdata))
    .translate_V_segment(V_alleles, intdata, V_segment)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### V_allele_has_stop_codon()
###

### Only needs access to the "coding_frame_start" column in 'intdata'.
### Returns a named logical vector that is parallel to 'V_alleles' and has
### the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'intdata' or for which 'intdata$coding_frame_start' has an NA.
V_allele_has_stop_codon <- function(V_alleles, intdata)
{
    V_aa <- translate_V_alleles(V_alleles, intdata)
    ans <- setNames(grepl("*", V_aa, fixed=TRUE), names(V_aa))
    ans[is.na(V_aa)] <- NA
    ans
}


### =========================================================================
### translate_codons() and related
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### remove_gaps()
###

GAP_LETTER <- "."

remove_gaps <- function(sequences, gap_letter=".")
{
    if (!(isSingleString(gap_letter) && nchar(gap_letter) == 1L))
        stop(wmsg("'gap_letter' must be a single letter"))
    if (is(sequences, "XStringSet")) {
        midx <- vmatchPattern(gap_letter, sequences, fixed=TRUE)
        at <- unname(as(midx, "CompressedIRangesList"))
    } else if (is(sequences, "XString")) {
        v <- matchPattern(gap_letter, sequences, fixed=TRUE)
        at <- as(v, "IRanges")
    } else {
        stop(wmsg("'sequences' must be an XStringSet or XString object"))
    }
    ## replaceAt() propagates the names and metadata on 'sequences'.
    replaceAt(sequences, at, value="")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_codons()
###

### Normalize a single offset.
.normarg_offset1 <- function(offset, dna_len)
{
    if (!isSingleNumber(offset))
        stop(wmsg("'offset' must be a single number"))
    if (!is.integer(offset))
        offset <- as.integer(offset)
    if (offset < 0L || offset > dna_len)
        stop(wmsg("'offset' must be >= 0 and <= 'length(dna)'"))
    offset
}

### Normalize a vector of offsets.
.normarg_offset2 <- function(offset, widths)
{
    if (!is.numeric(offset))
        stop(wmsg("'offset' must be numeric"))
    if (!is.integer(offset))
        offset <- as.integer(offset)
    noffset <- length(offset)
    nseq <- length(widths)
    if (noffset == 1L) {
        if (is.na(offset) || offset < 0L)
            stop(wmsg("'offset' cannot be NA or negative"))
        offset <- rep.int(offset, nseq)
        if (any(offset > widths))
            stop(wmsg("'offset' must be less than the length ",
                      "of the shortest sequence in 'dna'"))
    } else {
        if (noffset != nseq)
            stop(wmsg("'offset' must be a single value or ",
                      "a vector with one value per sequence in 'dna'"))
        if (anyNA(offset) || any(offset < 0L))
            stop(wmsg("'offset' cannot contain NAs or negative values"))
        if (any(offset > widths))
            stop(wmsg("each value in 'offset' must be less than the ",
                      "length of the corresponding sequence in 'dna'"))
    }
    offset
}

### Propagates names and metadata columns.
extract_codons <- function(dna, offset=0)
{
    if (is(dna, "DNAStringSet")) {
        dna_nchar <- nchar(dna)  # same as width(dna)
        offset <- .normarg_offset2(offset, dna_nchar)
    } else if (is(dna, "DNAString")) {
        dna_nchar <- nchar(dna)  # same as length(dna)
        offset <- .normarg_offset1(offset, dna_nchar)
    } else {
        stop(wmsg("'dna' must be a DNAStringSet or DNAString object"))
    }
    width <- dna_nchar - offset
    width <- width - width %% 3L
    subseq(dna, start=offset+1L, width=width)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_codons()
###

### Propagates names and metadata columns.
translate_codons <- function(dna, offset=0, with.init.codon=FALSE)
{
    if (!isTRUEorFALSE(with.init.codon))
        stop(wmsg("'with.init.codon' must be TRUE or FALSE"))
    codons <- extract_codons(dna, offset)
    ## Biostrings::translate() does not propagate the metadata columns at
    ## the moment but maybe it should?
    ans <- translate(codons, no.init.codon=!with.init.codon,
                             if.fuzzy.codon="solve")
    mcols(ans) <- mcols(codons)
    ans
}


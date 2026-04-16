### =========================================================================
### clean_allele_set() and related
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


checkarg_with.intdata <- function(with.intdata, gapped)
{
    if (!isTRUEorFALSE(with.intdata))
        stop(wmsg("'with.intdata' must be TRUE or FALSE"))
    if (with.intdata && !isTRUE(gapped))
        stop(wmsg("'with.intdata=TRUE' can only be used when 'gapped' is TRUE"))
}

.check_locus <- function(locus, what)
{
    if (!is.character(locus))
        stop(wmsg(what, " must be of type character"))
    if (anyNA(locus))
        stop(wmsg(what, " cannot contain NAs"))
    if (any(nchar(locus) != 3L))
        stop(wmsg(what, " must contain 3-letter strings"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Detection/handling of "repeated" vector elements
###
### By "repeated" vector elements here we mean elements with identical
### values **and** names. Note that, even though vector-like object 'x' will
### always be a DNAStringSet object in the context of our use cases, we've
### kept the implementation of these low-level utilities agnostic of the
### exact nature of 'x'.

.match_first_repeated <- function(x)
{
    x_names <- names(x)
    stopifnot(!is.null(x_names))
    m1 <- match(x, x)
    m2 <- match(x_names, x_names)
    selfmatchIntegerPairs(m1, m2)
}

### Similar to duplicated() but also takes into account the names on
### vector-like object 'x' to decide whether an element is "repeated" or not.
### Works on any vector-like object that has names and supports match().
.is_repeated <- function(x)
{
    .match_first_repeated(x) != seq_along(x)
}

### Alternate implementation.
### Note that:
### - This is a lot less generic than the implementation above because it
###   relies on splitAsList() working on 'x' and on duplicated() working on
###   the List derivative returned by splitAsList(). Works if 'x' is a
###   DNAStringSet object though, which is all we care about in the context
###   of this file.
### - Not as efficient as .is_repeated() above.
.is_repeated2 <- function(x)
{
    x_names <- names(x)
    stopifnot(!is.null(x_names))
    split_by_name <- splitAsList(x, x_names)
    unsplit(duplicated(split_by_name), x_names)
}

### MCOI: A character vector naming the Metadata Columns Of Interest.
.repeated_vector_elts_have_identical_metadata <- function(x, MCOI=NULL)
{
    stopifnot(is(x, "Vector"))
    x_mcols <- mcols(x, use.names=FALSE)
    stopifnot(is(x_mcols, "DataFrame"))
    if (!is.null(MCOI)) {
        stopifnot(is.character(MCOI))
        x_mcols <- x_mcols[ , MCOI, drop=FALSE]
    }
    m <- .match_first_repeated(x)
    identical(x_mcols[m, , drop=FALSE], x_mcols)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .drop_repeated_alleles()
###

### Has its own tests in tests/testthat/test-clean_allele_set.R!
### By "repeated" alleles here we mean alleles with identical sequences
### **and** names. Note that:
### - we keep alleles with identical sequences but different names;
### - we also keep alleles with identical names but different sequences.
.drop_repeated_alleles <- function(allele_set, verbose=FALSE)
{
    stopifnot(is(allele_set, "XStringSet") || is.character(allele_set))
    drop_idx <- which(.is_repeated(allele_set))
    drop_idx_len <- length(drop_idx)
    if (drop_idx_len == 0L)
        return(allele_set)  # no-op
    if (verbose) {
        in1string <- paste(names(allele_set)[drop_idx], collapse=", ")
        msg <- c("Deleted ", drop_idx_len, " \"repeated\" alleles from ",
                 "original set of ", length(allele_set), " alleles: ",
                 in1string, ".")
        note <- c("Note that alleles are considered \"repeated\" when ",
                  "they share the same ungapped sequence and name.")
        message("  o ", wmsg(msg, margin=4L), "\n",
                "    ", wmsg(note, margin=4L), "\n")
    }
    allele_set[-drop_idx]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .stop_on_ambiguous_allele_names()
###

.stop_on_ambiguous_allele_names <- function(allele_names)
{
    stopifnot(is.character(allele_names))
    ambiguous_allele_names <- unique(allele_names[duplicated(allele_names)])
    in1string <- paste(ambiguous_allele_names, collapse=", ")
    stop(wmsg("The following allele names are ambiguous: ", in1string, "."),
         "\n  ",
         wmsg("Use 'disambiguate.allele.names=TRUE' to disambiguate ",
              "them. (Also using 'verbose=TRUE' will display the ",
              "disambiguation details.)"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .disambiguate_allele_names()
###

### Similar to base::make.unique() but mangles with suffixes made of
### lowercase letters.
.make_pool_of_suffixes <- function(min_pool_size)
{
    max_pool_size <- (length(letters)**8 - 1) / (length(letters) - 1) - 1
    if (min_pool_size > max_pool_size)
        stop(wmsg("too many duplicate seq ids"))
    ans <- character(0)
    for (i in 1:7) {
        ans <- c(ans, mkAllStrings(letters, i))
        if (length(ans) >= min_pool_size)
            return(ans)
    }
    ## Should never happen because we checked for this condition earlier (see
    ## above).
    stop(wmsg("too many duplicate seq ids"))
}

.make_unique_seqids <- function(seqids)
{
    stopifnot(is.character(seqids))
    if (length(seqids) <= 1L)
        return(seqids)
    oo <- order(seqids)
    seqids2 <- seqids[oo]
    ir <- IRanges(1L, runLength(Rle(seqids2)))
    pool_of_suffixes <- .make_pool_of_suffixes(max(width(ir)))
    suffixes <- extractList(pool_of_suffixes, ir)  # CharacterList
    suffixes[lengths(suffixes) == 1L] <- ""
    seqids2 <- paste0(seqids2, unlist(suffixes, use.names=FALSE))
    ans <- seqids2[S4Vectors:::reverseIntegerInjection(oo, length(oo))]
    setNames(ans, names(seqids))
}

.disambiguate_allele_names <- function(allele_set, annotated=FALSE,
                                       verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))
    if (annotated) {
        allele_set_mcols <- mcols(allele_set, use.names=FALSE)
        stopifnot(is(allele_set_mcols, "DataFrame"),
                  identical(allele_names, allele_set_mcols[ , "allele_name"]))
    }
    if (anyDuplicated(allele_names)) {
        new_allele_names <- .make_unique_seqids(allele_names)
        if (verbose) {
            idx <- which(allele_names != new_allele_names)
            in1string <- paste0(allele_names[idx], "->", new_allele_names[idx],
                                collapse=", ")
            msg <- c("Disambiguation: renamed the ", length(idx), " ",
                     "following allele(s): ", in1string)
            message("  o ", wmsg(msg, margin=4L), "\n")
        }
        names(allele_set) <- new_allele_names
        if (annotated) {
            allele_set_mcols[ , "allele_name"] <- names(allele_set)
            mcols(allele_set, use.names=FALSE) <- allele_set_mcols
        }
    }
    allele_set
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .finish_cleanup()
###

.finish_cleanup <- function(allele_set, annotated=FALSE,
                            disambiguate.allele.names=FALSE,
                            verbose=FALSE)
{
    ## (B) Drop "repeated" alleles.
    allele_set <- .drop_repeated_alleles(allele_set, verbose=verbose)

    allele_names <- names(allele_set)
    if (anyDuplicated(allele_names)) {
        ## Note that, because we did (B), allele with repeated names
        ## are now guaranteed to have distinct DNA sequences.
        if (!disambiguate.allele.names)
            .stop_on_ambiguous_allele_names(allele_names)
        ## (C) Mangle allele names to make them unique if they're not.
        allele_set <- .disambiguate_allele_names(allele_set,
                                                 annotated=annotated,
                                                 verbose=verbose)
    }

    allele_set
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_V_allele_set() helpers
###

### 'allele_set' must have a "locus" metadata column.
### Returns a DNAStringSet object that:
### - is parallel to input DNAStringSet object 'allele_set';
### - carries the names and metadata columns of input object 'allele_set';
### - also carries the annotations obtained with compute_V_gene_delineations()
###   in additional metadata columns.
.annotate_V_alleles <- function(allele_set)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))

    allele_set_mcols <- mcols(allele_set, use.names=FALSE)
    if (is.null(allele_set_mcols) || is.null(locus <- allele_set_mcols$locus))
        stop(wmsg("'allele_set' must have a \"locus\" metadata column ",
                  "when 'with.intdata' is TRUE"))
    .check_locus(locus, "the \"locus\" metadata column on 'allele_set'")

    ## Annotate the V alleles based on their gaps.
    intdata <- compute_V_gene_delineations(allele_set)
    stopifnot(identical(intdata[ , "allele_name"], allele_names))

    ## Add "chain_type" column.
    intdata$chain_type <- make_chain_type("V", locus)

    mcols(allele_set) <- cbind(allele_set_mcols, intdata)
    allele_set
}

### Check that "repeated" V alleles (i.e. V alleles with identical
### **ungapped** DNA sequences **and** names) are annotated identically.
.stop_if_repeated_V_alleles_have_discordant_annotations <- function(allele_set)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    ## We only care about the annotations found in the Metadata Columns
    ## Of Interest listed in 'MCOI' when comparing the V allele annotations.
    MCOI <- setdiff(names(NDM_DATA_COL2CLASS), "chain_type")
    ok <- .repeated_vector_elts_have_identical_metadata(allele_set, MCOI=MCOI)
    if (!ok) {
        msg <- c("V alleles with identical ungapped sequences and names ",
                 "must have gaps that result in identical annotations (i.e. ",
                 "in identical FWR/CDR boundaries and \"coding_frame_start\")")
        stop(wmsg(msg))
    }
}

.stop_if_allele_sequences_have_gaps <- function(ngaps)
{
    stopifnot(is.integer(ngaps))
    bad_idx <- which(ngaps != 0L)
    if (length(bad_idx) == 0L)
        return()
    allele_names <- names(ngaps)
    stopifnot(!is.null(allele_names))
    first_bad_allele <- allele_names[[bad_idx[[1L]]]]
    stop(wmsg("Some V allele sequences have gaps (", length(bad_idx), " ",
              "out of ", length(ngaps), " input V allele sequences, e.g. ",
              "allele ", first_bad_allele, "). Use 'gapped=TRUE' if the ",
              "supplied V allele sequences are gapped."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_J_allele_set() helpers
###

### Propagate the names.
.extract_codon_starts <- function(headers)
{
    parsed_headers <- parse_imgt_fasta_headers(headers)
    codon_starts <- parsed_headers[ , "codon_start"]
    setNames(as.integer(codon_starts), names(codon_starts))
}

.annotate_J_alleles <- function(allele_set, codon_starts=NULL)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))

    allele_set_mcols <- mcols(allele_set, use.names=FALSE)
    if (is.null(allele_set_mcols) || is.null(locus <- allele_set_mcols$locus))
        stop(wmsg("'allele_set' must have a \"locus\" metadata column ",
                  "when 'with.auxdata' is TRUE"))
    .check_locus(locus, "the \"locus\" metadata column on 'allele_set'")

    ## Annotate the J alleles.
    auxdata <- compute_auxdata(allele_set, codon_starts=codon_starts,
                                           no.warning=TRUE)

    ## Sanity checks.
    stopifnot(identical(auxdata[ , "allele_name"], allele_names))
    expected_chain_type <- make_chain_type("J", locus)
    stopifnot(identical(auxdata$chain_type, expected_chain_type))

    mcols(allele_set) <- cbind(allele_set_mcols, auxdata)
    allele_set
}

### Check that "repeated" J alleles (i.e. J alleles with identical DNA
### sequences **and** names) are annotated identically.
.stop_if_repeated_J_alleles_have_discordant_annotations <- function(allele_set)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    ## We only care about the annotations found in the Metadata Columns
    ## Of Interest listed in 'MCOI' when comparing the V allele annotations.
    MCOI <- setdiff(names(AUXDATA_COL2CLASS), "chain_type")
    ok <- .repeated_vector_elts_have_identical_metadata(allele_set, MCOI=MCOI)
    if (!ok) {
        msg <- c("J alleles with identical sequences and names must have ",
                 "identical annotations")
        stop(wmsg(msg))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_allele_set()
### clean_V_allele_set()
### clean_J_allele_set()
###
### These are the workhorses behind create_region_db(), create_V_region_db(),
### and create_J_region_db().
###
### They all perform an **enhanced** version of the cleaning/editing
### implemented in the edit_imgt_file.pl script included in IgBLAST.
### More precisely, they all do the following:
###   (A) clean FASTA headers, as with original edit_imgt_file.pl script;
###   (B) drops repeated alleles;
###   (C) disambiguates allele names.
###
### In addition to (A) + (B) + (C) above, clean_V_allele_set() performs the
### following steps:
###   (aV) annotate the V alleles if 'gapped' and 'with.intdata' are TRUE;
###    (G) remove gaps from sequences if 'gapped' is TRUE, as with original
###        edit_imgt_file.pl script.
###
### In addition to (A) + (B) + (C) above, clean_J_allele_set() performs the
### following step:
###   (aJ) annotate the J alleles if 'with.auxdata' is TRUE.

### Generic cleaning. Use to clean allele sets of type D or C.
### Note that "repeated" alleles (i.e. alleles with identical sequences
### DNA sequences **and** names) are dropped.
### If, after dropping the "repeated" alleles, the names of the
### remaining alleles are not unique, then an error will be raised,
### unless 'disambiguate.allele.names' is set to TRUE, in which case
### the ambiguous allele names will be disambiguated.
clean_allele_set <- function(allele_set,
                             disambiguate.allele.names=FALSE,
                             verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers))
    stopifnot(isTRUEorFALSE(disambiguate.allele.names))
    stopifnot(isTRUEorFALSE(verbose))

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- clean_imgt_fasta_headers(headers)

    ## (B) + (C).
    .finish_cleanup(allele_set,
                    disambiguate.allele.names=disambiguate.allele.names,
                    verbose=verbose)
}

### Set 'gapped' to TRUE if the sequences in 'allele_set' are gapped (note
### that only V allele sequences are allowed to have gaps).
### Set 'with.intdata' to TRUE if the sequences in 'allele_set' are
### gapped **and** the associated internal data should be computed. Note
### that in this case 'allele_set' **must** carry a "locus" metadata column.
### Returns a DNAStringSet object that contains the ungapped allele sequences
### and carries the metadata columns of input object 'allele_set', if any.
### Furthermore, if 'with.intdata' is TRUE, the returned object will also
### carry the annotations obtained with compute_V_gene_delineations() in
### additional metadata columns.
clean_V_allele_set <- function(allele_set, gapped=FALSE, with.intdata=FALSE,
                               disambiguate.allele.names=FALSE,
                               verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers))
    stopifnot(isTRUEorFALSE(gapped))
    checkarg_with.intdata(with.intdata, gapped)
    stopifnot(isTRUEorFALSE(disambiguate.allele.names))
    stopifnot(isTRUEorFALSE(verbose))

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- clean_imgt_fasta_headers(headers)

    ngaps <- setNames(vcountPattern(GAP_LETTER, allele_set), names(allele_set))
    if (gapped) {
        if (with.intdata) {
            ## (aV) Annotate the V alleles.
            allele_set <- .annotate_V_alleles(allele_set)
        } else {
            warn_if_allele_sequences_have_no_gaps(ngaps)
        }
        allele_set <- remove_gaps(allele_set)  # (G)
        if (with.intdata)
            .stop_if_repeated_V_alleles_have_discordant_annotations(allele_set)
    } else {
        .stop_if_allele_sequences_have_gaps(ngaps)
    }

    ## (B) + (C).
    ### Note that V alleles are considered "repeated" if the have
    ### identical **ungapped** sequences DNA sequences and names.
    .finish_cleanup(allele_set, annotated=with.intdata,
                    disambiguate.allele.names=disambiguate.allele.names,
                    verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_J_allele_set()
###

clean_J_allele_set <- function(allele_set,
                               excluded_J_alleles=NULL,
                               with.auxdata=FALSE, imgt.fasta=FALSE,
                               disambiguate.allele.names=FALSE, verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers))
    stopifnot(isTRUEorFALSE(with.auxdata))
    stopifnot(isTRUEorFALSE(imgt.fasta))
    stopifnot(isTRUEorFALSE(disambiguate.allele.names))
    stopifnot(isTRUEorFALSE(verbose))

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- names(headers) <- clean_imgt_fasta_headers(headers)

    ## Drop excluded J alleles.
    if (!is.null(excluded_J_alleles)) {
        stopifnot(is.character(excluded_J_alleles))
        drop_idx <- which(names(allele_set) %in% excluded_J_alleles)
        if (length(drop_idx) != 0L) {
            allele_set <- allele_set[-drop_idx]
            headers <- headers[-drop_idx]
        }
    }

    if (with.auxdata) {
        if (imgt.fasta) {
            codon_starts <- .extract_codon_starts(headers)
        } else {
            codon_starts <- NULL
        }
        ## (aJ) Annotate the J alleles.
        allele_set <- .annotate_J_alleles(allele_set,
                                          codon_starts=codon_starts)
        .stop_if_repeated_J_alleles_have_discordant_annotations(allele_set)
    }

    ## (B) + (C).
    .finish_cleanup(allele_set, annotated=with.auxdata,
                    disambiguate.allele.names=disambiguate.allele.names,
                    verbose=verbose)
}


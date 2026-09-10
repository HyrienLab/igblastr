### =========================================================================
### clean_allele_set() and family
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


checkarg_auto.intdata <- function(auto.intdata, gapped)
{
    if (!isTRUEorFALSE(auto.intdata))
        stop(wmsg("'auto.intdata' must be TRUE or FALSE"))
    if (auto.intdata && !isTRUE(gapped))
        stop(wmsg("'auto.intdata=TRUE' can only be used when 'gapped' is TRUE"))
}

check_locus <- function(locus, what)
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
### .make_DORsummary()
###

.stop_on_ambiguous_allele_names <- function(allele_names)
{
    stopifnot(is.character(allele_names))
    is_dup <- duplicated(allele_names)
    ambiguous_allele_names <- unique(allele_names[is_dup])
    in1string <- paste(ambiguous_allele_names, collapse=", ")
    stop(wmsg("The following allele names are ambiguous: ", in1string, "."),
         "\n  ",
         wmsg("Use 'disambiguate.allele.names=TRUE' to disambiguate ",
              "them. (Also using 'verbose=TRUE' will display the ",
              "disambiguation details.)"))
}

### Has its own unit tests in tests/testthat/test-clean_allele_set.R!
### Make the "Drop Or Rename summary", a data.frame with one row per allele
### in 'allele_set', and the following columns:
###   - allele_name: the names on 'allele_set'.
###   - suffix: the disambiguation suffix (possible the empty string)
###     or NA if the allele should be dropped. The disambiguation suffix is
###     the suffix to add to the allele name in order to make allele names
###     unique across the allele set.
###   - locus: [optional] the "locus" metadata column from 'allele_set' if
###     the latter is an XStringSet derivative with such metadata column.
.make_DORsummary <- function(allele_set, disambiguate.allele.names=FALSE)
{
    ## The unit tests actually test .make_DORsummary() on character vectors!
    stopifnot(is(allele_set, "XStringSet") || is.character(allele_set))
    suffix <- rep.int(NA_character_, length(allele_set))
    remaining_idx <- which(!.is_repeated(allele_set))
    remaining_allele_names <- names(allele_set)[remaining_idx]
    if (anyDuplicated(remaining_allele_names) && !disambiguate.allele.names)
        .stop_on_ambiguous_allele_names(remaining_allele_names)
    suffix[remaining_idx] <- make_unique_allele_names(remaining_allele_names,
                                                      suffixes.only=TRUE)
    DORsummary <- data.frame(allele_name=names(allele_set), suffix=suffix)
    if (is(allele_set, "XStringSet"))
        DORsummary$locus <- mcols(allele_set)$locus
    DORsummary
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .drop_repeated_alleles()
###

### Has its own unit tests in tests/testthat/test-clean_allele_set.R!
### By "repeated" alleles we mean alleles with identical sequences **and**
### names. Note that:
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
        msg <- c("Dropped ", drop_idx_len, " \"repeated\" allele(s) from ",
                 "set of ", length(allele_set), " alleles: ", in1string, ".")
        note <- c("Note that alleles are considered \"repeated\" when ",
                  "they share the same ungapped sequence and name.")
        message("  o ", wmsg(msg, margin=4L), "\n",
                "    ", wmsg(note, margin=4L), "\n")
    }
    allele_set[-drop_idx]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .disambiguate_allele_names()
###

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
        new_allele_names <- make_unique_allele_names(allele_names)
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
### .drop_or_rename_alleles()
###

### Returns the modified 'allele_set' which is not necessarily parallel to
### the input 'allele_set' because some alleles may have been dropped.
### Also the returned object has the "Drop Or Rename summary" set on it (in
### its metadata() component). This is a 2-col data.frame that is **parallel**
### to the input 'allele_set'. See function .make_DORsummary() above in this
### file for the details.
.drop_or_rename_alleles <- function(allele_set, annotated=FALSE,
                                    disambiguate.allele.names=FALSE,
                                    verbose=FALSE)
{
    stopifnot(is(allele_set, "XStringSet"))
    DORsummary <- .make_DORsummary(allele_set,
                        disambiguate.allele.names=disambiguate.allele.names)

    ## (C) Drop "repeated" alleles.
    cleaned_allele_set <- .drop_repeated_alleles(allele_set, verbose=verbose)

    ## (D) Disambiguate allele names.
    allele_names <- names(cleaned_allele_set)
    if (anyDuplicated(allele_names)) {
        ## Note that .make_DORsummary() should already have caught this and
        ## raised an error if 'disambiguate.allele.names' is FALSE, but we
        ## check here anyways just in case one day .make_DORsummary() changes
        ## or we get rid of it.
        ## Also note that, because we did (C), alleles with repeated names
        ## are now guaranteed to have distinct DNA sequences.
        if (!disambiguate.allele.names)
            .stop_on_ambiguous_allele_names(allele_names)
        cleaned_allele_set <- .disambiguate_allele_names(cleaned_allele_set,
                                                         annotated=annotated,
                                                         verbose=verbose)
    }

    stopifnot(is.null(metadata(cleaned_allele_set)$DORsummary))
    metadata(cleaned_allele_set)$DORsummary <- DORsummary
    cleaned_allele_set
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
.annotate_V_alleles <-
    function(allele_set, fwrcdr_widths=IMGT_DEFAULT_FWRCDR_WIDTHS,
             verbose=FALSE)
{
    if (verbose) {
        msg <- c("Computing the intdata (a.k.a. igblastr-generated ",
                 "intdata) for the V alleles in the db")
        message("  o ", wmsg(msg, margin=4L), " ... ", appendLF=FALSE)
    }

    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))

    allele_set_mcols <- mcols(allele_set, use.names=FALSE)
    if (is.null(allele_set_mcols) || is.null(locus <- allele_set_mcols$locus))
        stop(wmsg("'allele_set' must have a \"locus\" metadata column ",
                  "when 'auto.intdata' is TRUE"))
    check_locus(locus, "the \"locus\" metadata column on 'allele_set'")

    ## Annotate the V alleles based on their gaps.
    intdata <- compute_V_gene_delineations(allele_set,
                                           fwrcdr_widths=fwrcdr_widths)
    stopifnot(identical(intdata[ , "allele_name"], allele_names))

    ## Add "chain_type" column.
    intdata$chain_type <- make_chain_type("V", locus)

    mcols(allele_set) <- cbind(allele_set_mcols, intdata)

    if (verbose)
        message("ok.\n")

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
    parsed_headers <- parse_imgt_fasta_headers(headers)  # character matrix
    codon_starts <- parsed_headers[ , "codon_start"]
    setNames(as.integer(codon_starts), names(headers))
}

.auxdata_completeness <- function(auxdata, verbose=FALSE)
{
    unsolved_idx <- which(is.na(auxdata[ , "cdr3_end"]))
    complete <- length(unsolved_idx) == 0L
    if (verbose) {
        if (complete) {
            in1string <- "none!"
        } else {
            unsolved_J_alleles <- auxdata[unsolved_idx, "allele_name"]
            in1string <- paste0(paste(unsolved_J_alleles, collapse=", "), ".")
        }
        msg <- c("Unsolved J allele(s): ", in1string)
        message("         --> ", wmsg(msg, margin=13L))
    }
    complete
}

### A verbose wrapper to compute_auxdata().
.compute_auxdata_STEP1 <- function(J_alleles, codon_starts, verbose=FALSE)
{
    if (verbose) {
        msg <- "Use compute_auxdata() to produce initial auxdata"
        message("    - STEP1: ", wmsg(msg, margin=13L), " ... ", appendLF=FALSE)
    }
    auxdata <- compute_auxdata(J_alleles, codon_starts=codon_starts,
                                          no.warnings=TRUE)
    if (verbose)
        message("ok.")
    auxdata
}

### A verbose wrapper to complete_auxdata_with_ref().
.compute_auxdata_STEP2 <- function(auxdata, ref_auxdata, verbose=FALSE)
{
    if (verbose) {
        ## We say that we're using validate_and_complete_auxdata_with_ref()
        ## even though we're actually using complete_auxdata_with_ref().
        ## But note that the former is just the user-exposed version of the
        ## latter.
        msg <- c("Use validate_and_complete_auxdata_with_ref() to validate ",
                 "the initial auxdata and try to complete it")
        message("    - STEP2: ", wmsg(msg, margin=13L), " ... ", appendLF=FALSE)
    }
    if (is.null(ref_auxdata)) {
        msg <- "SKIPPED!"
    } else {
        auxdata <- complete_auxdata_with_ref(auxdata, ref_auxdata,
                                             "igblastr-generated",
                                             "IgBLAST-provided")
        msg <- "ok."
    }
    if (verbose)
        message(msg)
    auxdata
}

### A verbose wrapper to solve_cdr3_ends_using_fwr4_aa_comparisons().
.compute_auxdata_STEP3 <- function(auxdata, J_alleles, verbose=FALSE)
{
    if (verbose) {
        msg <- c("Use solve_cdr3_ends_using_fwr4_aa_comparisons() to try ",
                 "to complete the auxdata obtained at previous step")
        message("    - STEP3: ", wmsg(msg, margin=13L), " ... ", appendLF=FALSE)
    }
    if (anyNA(auxdata[ , "cdr3_end"])) {
        auxdata <-
            solve_cdr3_ends_using_fwr4_aa_comparisons(auxdata, J_alleles)
        msg <- "ok."
    } else {
        msg <- "NOT NEEDED!"
    }
    if (verbose)
        message(msg)
    auxdata
}

### A verbose wrapper to solve_cdr3_ends_using_fwr4_dna_comparisons().
.compute_auxdata_STEP4 <- function(auxdata, J_alleles, verbose=FALSE)
{
    if (verbose) {
        msg <- c("Use solve_cdr3_ends_using_fwr4_dna_comparisons() to try ",
                 "to complete the auxdata obtained at previous step")
        message("    - STEP4: ", wmsg(msg, margin=13L), " ... ", appendLF=FALSE)
    }
    if (anyNA(auxdata[ , "cdr3_end"])) {
        auxdata <-
            solve_cdr3_ends_using_fwr4_dna_comparisons(auxdata, J_alleles)
        msg <- "ok."
    } else {
        msg <- "NOT NEEDED!"
    }
    if (verbose)
        message(msg)
    auxdata
}

### 'ref_auxdata' (if provided) is used to "complete" the computed
### auxiliary data. Note that the completion process will also verify
### that the computed auxiliary data is concordant with 'ref_auxdata'.
.annotate_J_alleles <- function(allele_set, codon_starts=NULL, ref_auxdata=NULL,
                                verbose=FALSE)
{
    if (verbose) {
        msg <- c("Computing the auxdata (a.k.a. igblastr-generated ",
                 "auxdata) for the J alleles in the db:")
        message("  o ", wmsg(msg, margin=4L))
    }

    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))

    allele_set_mcols <- mcols(allele_set, use.names=FALSE)
    if (is.null(allele_set_mcols) || is.null(locus <- allele_set_mcols$locus))
        stop(wmsg("'allele_set' must have a \"locus\" metadata column ",
                  "when 'auto.auxdata' is TRUE"))
    check_locus(locus, "the \"locus\" metadata column on 'allele_set'")

    auxdata <- .compute_auxdata_STEP1(allele_set, codon_starts, verbose=verbose)

    ## Sanity checks.
    stopifnot(identical(auxdata[ , "allele_name"], allele_names))
    expected_chain_type <- make_chain_type("J", locus)
    stopifnot(identical(auxdata[ , "chain_type"], expected_chain_type))

    auxdata <- .compute_auxdata_STEP2(auxdata, ref_auxdata, verbose=verbose)
    complete <- .auxdata_completeness(auxdata, verbose=verbose)

    if (!complete) {
        auxdata <- .compute_auxdata_STEP3(auxdata, allele_set, verbose=verbose)
        complete <- .auxdata_completeness(auxdata, verbose=verbose)
    }

    if (!complete) {
        auxdata <- .compute_auxdata_STEP4(auxdata, allele_set, verbose=verbose)
        complete <- .auxdata_completeness(auxdata, verbose=verbose)
    }

    mcols(allele_set) <- cbind(allele_set_mcols, auxdata)

    if (verbose)
        message("")

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
###   (C) drop "repeated" alleles;
###   (D) disambiguate allele names.
###
### In addition to (A) + (C) + (D) above, clean_V_allele_set() performs the
### following steps:
###   (Bv1) annotate the V alleles if 'gapped' and 'auto.intdata' are TRUE;
###   (Bv2) remove gaps from sequences if 'gapped' is TRUE, as with original
###         edit_imgt_file.pl script.
###
### In addition to (A) + (C) + (D) above, clean_J_allele_set() performs the
### following step:
###   (Bj) annotate the J alleles if 'auto.auxdata' is TRUE.

### Generic cleaning. Performs: (A) -> (C) -> (D).
### Use it to clean a set of germline D or C gene alleles.
### If, after performing the (C) step (i.e. drop the "repeated" alleles),
### the names of the remaining alleles are not unique, then an error will
### be raised, unless 'disambiguate.allele.names' is set to TRUE, in which
### case the ambiguous allele names will be disambiguated.
clean_allele_set <- function(allele_set,
                             disambiguate.allele.names=FALSE,
                             verbose=FALSE)
{
    stopifnot(is(allele_set, "XStringSet") || is.character(allele_set))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers),
              isTRUEorFALSE(disambiguate.allele.names),
              isTRUEorFALSE(verbose))

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- clean_imgt_fasta_headers(headers)

    ## (C) + (D).
    .drop_or_rename_alleles(allele_set,
                    disambiguate.allele.names=disambiguate.allele.names,
                    verbose=verbose)
}

### Performs: (A) -> (Bv1) -> (Bv2) -> (C) -> (D).
### Set 'gapped' to TRUE if the sequences in 'allele_set' are gapped (note
### that only V allele sequences are allowed to have gaps).
### Set 'auto.intdata' to TRUE if the sequences in 'allele_set' are
### gapped **and** the associated internal data should be computed. Note
### that in this case 'allele_set' **must** carry a "locus" metadata column.
### Returns a DNAStringSet object that contains the ungapped allele sequences
### and carries the metadata columns of input object 'allele_set', if any.
### Furthermore, if 'auto.intdata' is TRUE, the returned object will also
### carry the annotations obtained with compute_V_gene_delineations() in
### additional metadata columns.
clean_V_allele_set <- function(allele_set,
                               gapped=FALSE, auto.intdata=FALSE,
                               fwrcdr_widths=IMGT_DEFAULT_FWRCDR_WIDTHS,
                               disambiguate.allele.names=FALSE,
                               verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers),
              isTRUEorFALSE(gapped),
              isTRUEorFALSE(disambiguate.allele.names),
              isTRUEorFALSE(verbose))
    checkarg_auto.intdata(auto.intdata, gapped)

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- clean_imgt_fasta_headers(headers)

    ngaps <- setNames(vcountPattern(GAP_LETTER, allele_set), names(allele_set))
    if (gapped) {
        if (auto.intdata) {
            ## (Bv1) Annotate the V alleles.
            allele_set <- .annotate_V_alleles(allele_set,
                                              fwrcdr_widths=fwrcdr_widths,
                                              verbose=verbose)
        } else {
            warn_if_allele_sequences_have_no_gaps(ngaps)
        }
        allele_set <- remove_gaps(allele_set)  # (Bv2)
        if (auto.intdata)
            .stop_if_repeated_V_alleles_have_discordant_annotations(allele_set)
    } else {
        .stop_if_allele_sequences_have_gaps(ngaps)
    }

    ## (C) + (D).
    ### Note that V alleles are considered "repeated" if the have
    ### identical **ungapped** sequences DNA sequences and names.
    .drop_or_rename_alleles(allele_set, annotated=auto.intdata,
                    disambiguate.allele.names=disambiguate.allele.names,
                    verbose=verbose)
}

### Performs: (A) -> (Bj) -> (C) -> (D).
clean_J_allele_set <- function(allele_set, imgt.fasta.headers=FALSE,
                               auto.auxdata=FALSE, ref_auxdata=NULL,
                               disambiguate.allele.names=FALSE,
                               verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers),
              isTRUEorFALSE(imgt.fasta.headers),
              isTRUEorFALSE(auto.auxdata),
              isTRUEorFALSE(disambiguate.allele.names),
              isTRUEorFALSE(verbose))

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- clean_imgt_fasta_headers(headers)

    if (auto.auxdata) {
        if (imgt.fasta.headers) {
            names(headers) <- names(allele_set)
            codon_starts <- .extract_codon_starts(headers)
        } else {
            codon_starts <- NULL
        }
        ## (Bj) Annotate the J alleles.
        allele_set <- .annotate_J_alleles(allele_set,
                                          codon_starts=codon_starts,
                                          ref_auxdata=ref_auxdata,
                                          verbose=verbose)
        .stop_if_repeated_J_alleles_have_discordant_annotations(allele_set)
    }

    ## (C) + (D).
    .drop_or_rename_alleles(allele_set, annotated=auto.auxdata,
                    disambiguate.allele.names=disambiguate.allele.names,
                    verbose=verbose)
}


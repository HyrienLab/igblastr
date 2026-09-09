### =========================================================================
### Compute IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### For alleles in the IGHJ group (i.e. BCR germline J gene alleles on the
### heavy chain), the FWR4 region is expected to start with AA motif "WGXG.
.WGXG_motif <- "WG.G"
.WGXG_motif_dna <- DNAString("TGGGGNNNNGGN")  # reverse-translation of "WGXG"

### For all other J alleles, that is, for alleles in the IG[KL]J groups
### (i.e. BCR germline J gene alleles on the light chain) and all TCR
### germline J gene alleles, the FWR4 region is expected to start with
### AA motif "FGXG".
.FGXG_motif <- "FG.G"
.FGXG_motif_dna <- DNAString("TTYGGNNNNGGN")  # reverse-translation of "FGXG"


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .find_motif()
###

### Returns an integer vector parallel to 'J_alleles' that contains the
### 1-based start position of the motif for each sequence in 'J_alleles'.
### For alleles that don't contain the motif, the position is set to NA.
### For alleles that contain more than one occurrence of the motif, we
### return the position of the first occurrence.
.find_motif_dna <- function(J_dna, motif_dna)
{
    stopifnot(is(J_dna, "DNAStringSet"), is(motif_dna, "DNAString"))
    m <- vmatchPattern(motif_dna, J_dna, fixed=FALSE)
    as.integer(heads(start(m), n=1L))
}

.find_motif_aa <- function(J_aa, motif)
{
    stopifnot(is(J_aa, "AAStringSet"), isSingleNonWhiteString(motif))
    m <- gregexpr(motif, as.character(J_aa))
    ans <- as.integer(heads(m, n=1L))
    ans[ans == -1L] <- NA_integer_
    ans
}

.find_motif <- function(J_alleles, codon_starts, motif_dna, motif)
{
    stopifnot(is(J_alleles, "DNAStringSet"), is.integer(codon_starts),
              length(J_alleles) == length(codon_starts))
    ans <- rep.int(NA_integer_, length(J_alleles))
    idx1 <- which(is.na(codon_starts))
    if (length(idx1) != 0L)
        ans[idx1] <- .find_motif_dna(J_alleles[idx1], motif_dna)
    idx2 <- which(!is.na(codon_starts))
    if (length(idx2) != 0L) {
        codon_starts2 <- codon_starts[idx2]
        J_aa <- translate_codons(J_alleles[idx2], offset=codon_starts2 - 1L)
        starts2 <- .find_motif_aa(J_aa, motif)
        ans[idx2] <- codon_starts2 + 3L * (starts2 - 1L)
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_auxdata()
###

.normarg_codon_starts <- function(codon_starts, allele_names, fasta_headers)
{
    if (is.null(codon_starts)) {
        codon_starts <- setNames(rep.int(NA_integer_, length(allele_names)),
                                 allele_names)
        return(codon_starts)
    }
    if (!is.numeric(codon_starts))
        stop(wmsg("'codon_starts' must be NULL or an integer vector"))
    if (length(codon_starts) != length(allele_names))
        stop(wmsg("'codon_starts' must have the same length as 'J_alleles'"))
    codon_starts_names <- names(codon_starts)
    if (is.null(codon_starts_names))
        stop(wmsg("'codon_starts' must carry the names of the J alleles"))
    ## Because of clean_imgt_fasta_headers()'s funky business with the
    ## IG allele names for Mus spretus, the allele names in 'allele_names'
    ## are not necessarily the same as in 'fasta_headers[ , "allele_name"]'
    ## (some allele names in the former can have the _Mus_spretus suffix
    ## added to them), so we compare 'codon_starts_names' to both.
    if (!identical(codon_starts_names, allele_names)) {
        parsed_headers <- parse_imgt_fasta_headers(fasta_headers)
        if (!identical(codon_starts_names, parsed_headers[ , "allele_name"]))
            stop(wmsg("'codon_starts' must carry the names of the J alleles, ",
                      "and be in the same order as 'J_alleles'"))
        names(codon_starts) <- allele_names
    }
    if (!is.integer(codon_starts))
        codon_starts <- setNames(as.integer(codon_starts), names(codon_starts))
    if (!(all(codon_starts %in% c(1:3, NA_integer_))))
        stop(wmsg("the non-NA values in 'codon_starts' must ",
                  "be >= 1 and <= 3"))
    codon_starts
}

### Not exported!
J_alleles_with_missing_coding_frame_start <- function(auxdata)
{
    stopifnot(is.data.frame(auxdata))
    bad_idx <- which(is.na(auxdata[ , "coding_frame_start"]))
    auxdata[bad_idx, "allele_name"]
}

### Warn user if "coding frame start" could not be determined for some alleles.
.warn_if_missing_coding_frame_starts <- function(auxdata)
{
    bad_alleles <- J_alleles_with_missing_coding_frame_start(auxdata)
    if (length(bad_alleles) != 0L) {
        in1string <- paste(bad_alleles, collapse=", ")
        warning(wmsg("The \"coding frame start\" could not be determined ",
                     "for J allele(s): ", in1string, "."),
                "\n  ",
                wmsg("--> coding_frame_start, cdr3_end, and extra_bps ",
                     "were set to NA for these alleles"))
    }
}

### Not exported!
### Very rare but can happen (we've only seen this for rhesus monkey
### allele IGHJ0-ZXTW*01 from OGRDB so far).
warn_if_negative_cdr3_end <- function(auxdata, what)
{
    bad_idx <- which(auxdata[ , "cdr3_end"] < 0L)
    if (length(bad_idx) != 0L) {
        in1string <- paste(auxdata[bad_idx , "allele_name"], collapse=", ")
        warning(wmsg("The ", what, " is negative for J allele(s): ",
                     in1string, "."),
                "\n  ",
                wmsg("Note that you won't be able to save the returned ",
                     "data.frame with write_auxdata() unless you replace ",
                     "the negative \"cdr3_end\" value(s) with NA(s) or drop ",
                     "the row(s) with a negative \"cdr3_end\" value."))
    }
}

### Compute the auxiliary data for 'J_alleles' by searching motifs WGXG
### and FGXG.
### Returns a data.frame with 1 row per sequence in 'J_alleles' and with
### the same columns as the data.frame returned by load_auxdata() (see
### file R/load_auxdata.R).
compute_auxdata <- function(J_alleles, codon_starts=NULL, no.warnings=FALSE)
{
    if (!is(J_alleles, "DNAStringSet"))
        stop(wmsg("'J_alleles' must be a DNAStringSet object"))
    fasta_headers <- names(J_alleles)
    if (is.null(fasta_headers))
        stop(wmsg("'J_alleles' must have names"))

    allele_names <- clean_imgt_fasta_headers(fasta_headers)
    allele_loci <- substr(allele_names, 1L, 3L)
    loci_prefix <- extract_loci_prefix(allele_loci)
    if (is.na(loci_prefix))
        stop(wmsg("all allele names must start either ",
                  "with 'IG[HKL]' or with 'TR[ABGD]'"))
    if (!all(substr(allele_names, 4L, 4L) == "J"))
        stop(wmsg("the 4th letter in every J allele name must be \"J\""))

    codon_starts <- .normarg_codon_starts(codon_starts, allele_names,
                                          fasta_headers)

    if (!isTRUEorFALSE(no.warnings))
        stop(wmsg("'no.warnings' must be TRUE or FALSE"))

    ## Compute 'cdr3_ends'.
    ## We use WGXG motif for J alleles from IGH locus, and FGXG motif for
    ## J alleles from all other loci.
    fwr4_starts <- rep.int(NA_integer_, length(J_alleles))
    idx1 <- which(allele_loci == "IGH")
    fwr4_starts[idx1] <- .find_motif(J_alleles[idx1], codon_starts[idx1],
                                     .WGXG_motif_dna, .WGXG_motif)
    idx2 <- which(allele_loci != "IGH")
    fwr4_starts[idx2] <- .find_motif(J_alleles[idx2], codon_starts[idx2],
                                     .FGXG_motif_dna, .FGXG_motif)
    cdr3_ends <- fwr4_starts - 1L  # 1-based

    ## Compute 'coding_frame_starts'.
    coding_frame_starts <- unname(codon_starts) - 1L
    stopifnot(all(coding_frame_starts == cdr3_ends %% 3L, na.rm=TRUE))
    idx <- which(is.na(coding_frame_starts) & !is.na(cdr3_ends))
    coding_frame_starts[idx] <- cdr3_ends[idx] %% 3L

    ## Compute 'chain_types' and 'extra_bps'.
    chain_types <- make_chain_type("J", allele_loci)
    extra_bps <- (width(J_alleles) - coding_frame_starts) %% 3L

    ans <- data.frame(
        allele_name       =allele_names,
        coding_frame_start=coding_frame_starts,  # 0-based
        chain_type        =chain_types,
        cdr3_end          =cdr3_ends - 1L,  # 0-based
        extra_bps         =extra_bps
    )

    if (!no.warnings) {
        .warn_if_missing_coding_frame_starts(ans)
        what <- "computed CDR3 end position"
        warn_if_negative_cdr3_end(ans, what)
    }

    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_discordant_auxdata()
### complete_auxdata_with_ref()
### validate_and_complete_auxdata_with_ref()
###

### Not exported!
### Has its own unit tests in tests/testthat/test-compute_auxdata.R!
### See find_discordant_rows() in R/utils.R for what it returns.
find_discordant_auxdata <- function(auxdata, ref_auxdata)
{
    check_auxdata_col2class(auxdata)
    check_auxdata_col2class(ref_auxdata)
    find_discordant_rows(auxdata, ref_auxdata, "allele_name")
}

### Not exported!
### Has its own unit tests in tests/testthat/test-compute_auxdata.R!
### See complete_df_with_ref() in R/utils.R for what it returns.
complete_auxdata_with_ref <- function(auxdata, ref_auxdata,
                                      auxdata_label="computed",
                                      ref_label="reference")
{
    check_auxdata_col2class(auxdata)
    check_auxdata_col2class(ref_auxdata)
    complete_df_with_ref(auxdata, ref_auxdata, "allele_name",
                         "\"auxiliary data\"", auxdata_label, ref_label)
}

### User-exposed wrapper to complete_auxdata_with_ref().
validate_and_complete_auxdata_with_ref <- function(auxdata, ref_auxdata)
    complete_auxdata_with_ref(auxdata, ref_auxdata)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The .get_fwr4_start_from_fwr4refset_*() functions
###

.normarg_max.dist <- function(max.dist, max_max.dist)
{
    if (!isSingleNumber(max.dist))
        stop(wmsg("'max.dist' must be a single integer"))
    if (!is.integer(max.dist))
        max.dist <- as.integer(max.dist)
    if (!(max.dist >= 0L && max.dist <= max_max.dist))
        stop(wmsg("'max.dist' must be >= 0 and <= ", max_max.dist))
    max.dist
}

.normarg_standout.by <- function(standout.by, max_standout.by)
{
    if (!isSingleNumber(standout.by))
        stop(wmsg("'standout.by' must be a single integer"))
    if (!is.integer(standout.by))
        standout.by <- as.integer(standout.by)
    if (!(standout.by >= 1L && standout.by <= max_standout.by))
        stop(wmsg("'standout.by' must be >= 1 and <= ", max_standout.by))
    standout.by
}

.normarg_min.score <- function(min.score)
{
    if (!isSingleNumber(min.score))
        stop(wmsg("'min.score' must be a single number"))
    if (!(min.score >= 0 && min.score <= 1))
        stop(wmsg("'min.score' must be >= 0 and <= 1"))
    min.score
}

.normarg_standout.by2 <- function(standout.by)
{
    if (!isSingleNumber(standout.by))
        stop(wmsg("'standout.by' must be a single number"))
    if (!(standout.by > 0 && standout.by <= 1))
        stop(wmsg("'standout.by' must be > 0 and <= 1"))
    standout.by
}

### 'fwr4refset' is a set of known FWR4 sequences that we use as reference.
### The function moves a sliding window of width 'window_width' to the
### positions in 'window_positions' on 'J_allele', and, for each position
### of the window, computes the Hamming distance between the sequence in the
### window and 'fwr4refset'.
### Returns the distances in an integer vector parallel to 'window_positions'.
.compute_sliding_window_dists_to_fwr4refset <-
    function(J_allele, window_positions, window_width,
             fwr4refset, min_fwr4_width)
{
    stopifnot(is(J_allele, "XString"),
              is.integer(window_positions), isSingleInteger(window_width),
              max(window_positions) + window_width - 1L <= nchar(J_allele),
              is(fwr4refset, "XStringSet"), isSingleInteger(min_fwr4_width),
              all(width(fwr4refset) >= min_fwr4_width))
    vapply(window_positions,
        function(pos) {
            window_seq <- subseq(J_allele, start=pos, width=window_width)
            min(neditAt(window_seq, fwr4refset))
        },
        integer(1)
    )
}

.get_fwr4_start_from_fwr4refset_dist <-
    function(J_allele, window_positions, window_width,
             fwr4refset, min_fwr4_width,
             max.dist, standout.by)
{
    stopifnot(is(J_allele, "XString"),
              is.integer(window_positions), length(window_positions) >= 1L,
              isSingleInteger(window_width),
              max(window_positions) + window_width - 1L <= nchar(J_allele),
              is(fwr4refset, "XStringSet"), isSingleInteger(min_fwr4_width),
              all(width(fwr4refset) >= min_fwr4_width),
              isSingleInteger(max.dist), isSingleInteger(standout.by))
    dists <- .compute_sliding_window_dists_to_fwr4refset(J_allele,
                                              window_positions, window_width,
                                              fwr4refset, min_fwr4_width)
    i <- which.min(dists)
    DP <- dists[[i]]
    if (DP > max.dist)
        return(NA_integer_)  # no P achieves D(P) <= max.dist
    if (sum(dists - standout.by < DP) != 1L)
        return(NA_integer_)  # P does not stand out in the crowd
    window_positions[[i]]    # returns P
}

### All FWR4 we've seen so far have at least 10 amino acids except in
### mouse J allele TRAJ45*02 from IMGT where it has only 9 amino acids
### (this allele sequence is probably truncated).
.SLIDING_WINDOW_AA_WIDTH <- 10L
.MIN_FWR4_AA_WIDTH <- 9L
.AA_MAX_MAX_DIST <- 5L
.AA_MAX_STANDOUT_BY <- 10L

### Has its own unit tests in tests/testthat/test-compute_auxdata.R!
### The "best FWR4 start position" in 'aa_string' is defined as follow:
### If D(p) is the distance between the 10-amino acid window in 'aa_string'
### that starts at position 'p' and the set of known FWR4 sequences
### in 'fwr4refset', then the "best FWR4 start position" is the position 'P'
### in 'aa_string' that minimize D(P). With the additional constraints that:
###   (a) D(P) <= 2
###   (b) D(p) >= D(P) + 5 for any other position 'p'
### In other words, P must stand out in the crowd!
.get_fwr4_start_from_fwr4refset_aa_dist <-
    function(aa_string, fwr4refset, max.dist=2, standout.by=5)
{
    max.dist <- .normarg_max.dist(max.dist, .AA_MAX_MAX_DIST)
    standout.by <- .normarg_standout.by(standout.by, .AA_MAX_STANDOUT_BY)
    npositions <- nchar(aa_string) - .SLIDING_WINDOW_AA_WIDTH + 1L
    if (npositions <= 0L)
        return(NA_integer_)
    window_positions <- seq_len(npositions)
    .get_fwr4_start_from_fwr4refset_dist(
                    aa_string, window_positions, .SLIDING_WINDOW_AA_WIDTH,
                    fwr4refset, .MIN_FWR4_AA_WIDTH,
                    max.dist, standout.by)
}

.SLIDING_WINDOW_DNA_WIDTH <- .SLIDING_WINDOW_AA_WIDTH * 3L
.MIN_FWR4_DNA_WIDTH <- .MIN_FWR4_AA_WIDTH * 3L
.DNA_MAX_MAX_DIST <- .AA_MAX_MAX_DIST * 3L
.DNA_MAX_STANDOUT_BY <- .AA_MAX_STANDOUT_BY * 3L

### Has its own unit tests in tests/testthat/test-compute_auxdata.R!
### Moves the sliding window to positions 1, 4, 7, etc..  on 'dna_string'.
.get_fwr4_start_from_fwr4refset_dna_dist <-
    function(dna_string, fwr4refset, max.dist=4, standout.by=15)
{
    max.dist <- .normarg_max.dist(max.dist, .DNA_MAX_MAX_DIST)
    standout.by <- .normarg_standout.by(standout.by, .DNA_MAX_STANDOUT_BY)
    npositions <- (nchar(dna_string) - .SLIDING_WINDOW_DNA_WIDTH) %/% 3L + 1L
    if (npositions <= 0L)
        return(NA_integer_)
    window_positions <- seq_len(npositions) * 3L - 2L
    .get_fwr4_start_from_fwr4refset_dist(
                    dna_string, window_positions, .SLIDING_WINDOW_DNA_WIDTH,
                    fwr4refset, .MIN_FWR4_DNA_WIDTH,
                    max.dist, standout.by)
}

### Has its own unit tests in tests/testthat/test-compute_auxdata.R!
.get_fwr4_start_from_fwr4refset_dna_PWM <-
    function(dna_string, pwm, min.score=0.80, standout.by=0.40)
{
    stopifnot(is.matrix(pwm),
              identical(rownames(pwm), c("A", "C", "G", "T")),
              ncol(pwm) == .SLIDING_WINDOW_DNA_WIDTH)
    min.score <- .normarg_min.score(min.score)
    standout.by <- .normarg_standout.by2(standout.by)
    npositions <- (nchar(dna_string) - .SLIDING_WINDOW_DNA_WIDTH) %/% 3L + 1L
    if (npositions <= 0L)
        return(NA_integer_)
    window_positions <- seq_len(npositions) * 3L - 2L
    scores <- PWMscoreStartingAt(pwm, dna_string, starting.at=window_positions)
    i <- which.max(scores)
    SP <- scores[[i]]
    if (SP < min.score)
        return(NA_integer_)  # no P achieves S(P) >= min.score
    if (sum(scores + standout.by > SP) != 1L)
        return(NA_integer_)  # P does not stand out in the crowd
    window_positions[[i]]    # returns P
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The .compute_fwr4_starts_from_fwr4_*() functions
###

### Returns an integer vector parallel to 'aa_strings' that contains the
### 1-based FWR4 start positions for each amnino acid sequence.
.get_fwr4_starts_from_fwr4refset_aa_dist <-
    function(aa_strings, fwr4refset, max.dist=2, standout.by=5)
{
    stopifnot(is(aa_strings, "AAStringSet"),
              is(fwr4refset, "AAStringSet"),
              all(nchar(fwr4refset) >= .MIN_FWR4_AA_WIDTH))
    vapply(seq_along(aa_strings),
        function(i)
            .get_fwr4_start_from_fwr4refset_aa_dist(
                            aa_strings[[i]], fwr4refset,
                            max.dist=max.dist, standout.by=standout.by),
        integer(1)
    )
}

### All .compute_fwr4_starts_from_fwr4_*() functions return an integer vector
### parallel to 'J_alleles' that contains the 1-based FWR4 start positions.

.compute_fwr4_starts_from_fwr4_aa_comparisons <-
    function(J_alleles, fwr4refset, max.dist=2, standout.by=5)
{
    stopifnot(is(J_alleles, "DNAStringSet"),
              is(fwr4refset, "DNAStringSet"),
              all(nchar(fwr4refset) >= .MIN_FWR4_DNA_WIDTH))
    aa_strings <- translate_codons(J_alleles)
    fwr4refset <- translate_codons(fwr4refset)
    aa_starts <- .get_fwr4_starts_from_fwr4refset_aa_dist(
                                  aa_strings, fwr4refset,
                                  max.dist=max.dist,
                                  standout.by=standout.by)
    aa_starts * 3L - 2L
}

.compute_fwr4_starts_from_fwr4_dna_comparisons <-
    function(J_alleles, fwr4refset, max.dist=4, standout.by=15)
{
    stopifnot(is(J_alleles, "DNAStringSet"),
              is(fwr4refset, "DNAStringSet"),
              all(nchar(fwr4refset) >= .MIN_FWR4_DNA_WIDTH))
    vapply(seq_along(J_alleles),
        function(i)
            .get_fwr4_start_from_fwr4refset_dna_dist(
                            J_alleles[[i]], fwr4refset,
                            max.dist=max.dist, standout.by=standout.by),
        integer(1)
    )
}

.build_fwr4_dna_PWM <- function(fwr4refset)
{
    stopifnot(is(fwr4refset, "DNAStringSet"),
              all(nchar(fwr4refset) >= .MIN_FWR4_DNA_WIDTH))
    ## We want to use a sliding window of width 30 so we need a PWM of
    ## width 30. This means that all the FWR4 sequences we use to build
    ## the PWM must have at least 30 nucleotides.
    ## The only FWR4 we know that is shorter than that is mouse J allele
    ## TRAJ45*02 from IMGT, which has only 27 nucleotides:
    ##    TTTGGGAAAGGAACTCAGCTGATCATC         # 27 nucleotides
    ## This seems to be a truncated version of the FWR4 of mouse J allele
    ## TRAJ45*01, which is:
    ##    TTTGGGAAAGGAACTCAGCTGATCATCCAGCCCT  # 34 nucleotides
    ## So we replace the former with the latter.
    bad_idx <- which(width(fwr4refset) < .SLIDING_WINDOW_DNA_WIDTH)
    if (length(bad_idx) != 0L) {
        ## First 30 nucleotides only.
        mouse_TRAJ45star01_fwr4 <- DNAString("TTTGGGAAAGGAACTCAGCTGATCATCCAG")
        stopifnot(nchar(mouse_TRAJ45star01_fwr4) == .SLIDING_WINDOW_DNA_WIDTH)
        ## We check that the short FWR4 sequences are prefixes
        ## of 'mouse_TRAJ45star01_fwr4'.
        prefixes <- Views(mouse_TRAJ45star01_fwr4,
                          start=1L,
                          end=width(fwr4refset)[bad_idx])
        stopifnot(all(fwr4refset[bad_idx] == as(prefixes, "DNAStringSet")))
        fwr4refset[bad_idx] <- DNAStringSet(mouse_TRAJ45star01_fwr4)
    }
    PWM(subseq(fwr4refset, start=1L, end=.SLIDING_WINDOW_DNA_WIDTH))
}

.compute_fwr4_starts_from_fwr4_dna_PWM <-
    function(J_alleles, fwr4refset, min.score=0.80, standout.by=0.40)
{
    stopifnot(is(J_alleles, "DNAStringSet"),
              is(fwr4refset, "DNAStringSet"),
              all(nchar(fwr4refset) >= .MIN_FWR4_DNA_WIDTH))
    pwm <- .build_fwr4_dna_PWM(fwr4refset)
    vapply(seq_along(J_alleles),
        function(i)
            .get_fwr4_start_from_fwr4refset_dna_PWM(
                            J_alleles[[i]], pwm,
                            min.score=min.score, standout.by=standout.by),
        integer(1)
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The solve_cdr3_ends_using_fwr4_*() functions
###

.solve_cdr3_ends <- function(auxdata, J_alleles, FUN, ...)
{
    check_auxdata_col2class(auxdata)
    if (!is(J_alleles, "DNAStringSet") || length(J_alleles) != nrow(auxdata))
        stop(wmsg("'J_alleles' must be a DNAStringSet object ",
                  "with one sequence per row in 'auxdata'"))
    fasta_headers <- names(J_alleles)
    if (is.null(fasta_headers))
        stop(wmsg("'J_alleles' must have names"))
    names(J_alleles) <- clean_imgt_fasta_headers(fasta_headers)
    if (!identical(names(J_alleles), auxdata[ , "allele_name"]))
        stop(wmsg("the allele names on 'J_alleles' must ",
                  "be identical to 'auxdata$allele_name'"))
    coding_frame_starts <- auxdata[ , "coding_frame_start"]
    if (anyNA(coding_frame_starts))
        stop(wmsg("'auxdata$coding_frame_start' cannot contain NAs"))

    cdr3_ends <- auxdata[ , "cdr3_end"]  # 0-based
    solve_me <- is.na(cdr3_ends)
    if (!any(solve_me))
        return(auxdata)

    unsolved_offsets <- coding_frame_starts[solve_me]
    unsolved_J_alleles <- subseq(J_alleles[solve_me],
                                 start=unsolved_offsets + 1L)
    fwr4refset <- subseq(J_alleles[!solve_me],
                         start=cdr3_ends[!solve_me] + 2L)
    solved_fwr4_starts <- FUN(unsolved_J_alleles, fwr4refset, ...)
    solved_cdr3_ends <- solved_fwr4_starts + unsolved_offsets - 2L  # 0-based
    auxdata[solve_me, "cdr3_end"] <- solved_cdr3_ends
    warn_if_negative_cdr3_end(auxdata, "CDR3 end position")
    auxdata
}

solve_cdr3_ends_using_fwr4_aa_comparisons <-
    function(auxdata, J_alleles, max.dist=2, standout.by=5)
{
    .solve_cdr3_ends(auxdata, J_alleles,
                     .compute_fwr4_starts_from_fwr4_aa_comparisons,
                     max.dist=max.dist, standout.by=standout.by)
}

solve_cdr3_ends_using_fwr4_dna_comparisons <-
    function(auxdata, J_alleles, max.dist=4, standout.by=15)
{
    .solve_cdr3_ends(auxdata, J_alleles,
                     .compute_fwr4_starts_from_fwr4_dna_comparisons,
                     max.dist=max.dist, standout.by=standout.by)
}

solve_cdr3_ends_using_fwr4_dna_PWM <-
    function(auxdata, J_alleles, min.score=0.80, standout.by=0.40)
{
    .solve_cdr3_ends(auxdata, J_alleles,
                     .compute_fwr4_starts_from_fwr4_dna_PWM,
                     min.score=min.score, standout.by=standout.by)
}


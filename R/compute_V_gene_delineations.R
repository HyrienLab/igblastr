### =========================================================================
### compute_V_gene_delineations()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FWR/CDR fixed end positions w.r.t. the V gapped sequences
###

### The IMGT unique numbering provides a standardized delimitation of
### the FWR and CDR regions. This standard relies on the maximum observed
### FWR/CDR widths reported by IMGT. Gaps are inserted in the germline V gene
### protein sequences so that the widths of the gapped FWR/CDR regions are
### effectively the maximum observed widths.
### See https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
IMGT_DEFAULT_FWRCDR_WIDTHS <-
    c(fwr1=26L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L)

### There are some exceptions to this rule though.

### Not exported!
IMGT_MOUSE_FWRCDR_WIDTHS <- rbind(
    ## Not 100% sure about the CDR2/FWR3 junction for the IMGT V alleles
    ## on loci TRA and TRD, hence our hesitation between using
    ## cdr2=11L/fwr3=39L and cdr2=10L/fwr3=40L. However, using the latter
    ## is A LOT MORE in agreement with 'load_intdata("mouse")'.
    ## More precisely, with the former, doing
    ##   install_IMGT_germline_db("202614-2", "Mus_musculus", tcr.db=TRUE)
    ## introduces CDR2/FWR3 junction disagreements for 150+ TRA alleles
    ## and 8 TRD alleles. While with the latter, we get CDR2/FWR3 junction
    ## disagreements for 0 TRA allele and only 3 TRD alleles: TRDV5*01,
    ## TRDV5*03, and TRDV5*04!
    #TRA=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=11L, fwr3=39L),
    TRA=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L),
    #TRD=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=11L, fwr3=39L),
    TRD=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L)
)

### Not exported!
IMGT_RHESUS_MONKEY_FWRCDR_WIDTHS <- rbind(
    IGH=c(fwr1=27L, cdr1=13L, fwr2=17L, cdr2=10L, fwr3=39L),
    IGK=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L),
    IGL=c(fwr1=27L, cdr1=12L, fwr2=19L, cdr2=10L, fwr3=39L),
    ## Not sure about the CDR2/FWR3 junction for the IMGT V alleles
    ## on loci TRA, TRB, and TRG, hence our hesitation between using
    ## cdr2=11L/fwr3=39L and cdr2=10L/fwr3=40L, even though the latter
    ## seems more likely to be "the truth" (gut feeling based on some
    ## observations that are too long to explain here). Note that we cannot
    ## disambiguate by comparing with 'load_intdata("rhesus_monkey")' like
    ## we did for IMGT_MOUSE_FWRCDR_WIDTHS above because IgBLAST does not
    ## provide internal data for rhesus monkey TR alleles. So for now, we
    ## disable automatic intdata generation in
    ##   install_IMGT_germline_db("<release>", "Macaca_mulatta", tcr.db=TRUE)
    ## See .from_auto.intdata_to_intdata() in R/install_IMGT_germline_db.R.
    #TRA=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=11L, fwr3=39L),
    TRA=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L),
    #TRB=c(fwr1=26L, cdr1=12L, fwr2=17L, cdr2=11L, fwr3=39L),
    TRB=c(fwr1=26L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L),
    #TRG=c(fwr1=26L, cdr1=12L, fwr2=17L, cdr2=11L, fwr3=39L),
    TRG=c(fwr1=26L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L)
)

### Not exported!
### In the IMGT gapped sequences for rainbow trout, the end of the FWR1 regions
### is at position 29 (in amino acid space) instead of standard position 26. As
### a result, the CDR1/FWR2/CDR2/FWR3 are shifted downstream by 3 positions.
IMGT_RAINBOW_TROUT_FWRCDR_WIDTHS <-
    c(fwr1=29L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normarg_gapped_V_alleles()
### normarg_fwrcdr_widths()
### warn_if_allele_sequences_have_no_gaps()
###

### Not exported!
### Input must be the path to a FASTA file or a named DNAStringSet object.
### If the former, returns a BStringSet object. Otherwise, returns the
### input DNAStringSet object. In both cases, the names on the returned
### object are passed thru clean_imgt_fasta_headers().
normarg_gapped_V_alleles <- function(gapped_V_alleles)
{
    if (isSingleString(gapped_V_alleles)) {
        what <- paste0("some allele names in ", gapped_V_alleles)
        ## Some IMGT FASTA files (e.g. for Aotus_nancymaae and
        ## Nonhuman_primates) have nucleotide sequences that contain
        ## the letter 'x'. Not sure what that's supposed to represent.
        ## Note that a well established consensus is to use 'n' or 'N' to
        ## represent an unknown nucleotide (wildcard). Anyways, this breaks
        ## readDNAStringSet() so we use readBStringSet() instead.
        gapped_V_alleles <- readBStringSet(gapped_V_alleles)
        allele_names <- names(gapped_V_alleles)
    } else if (is(gapped_V_alleles, "DNAStringSet")) {
        allele_names <- names(gapped_V_alleles)
        if (is.null(allele_names))
            stop(wmsg("DNAStringSet object 'gapped_V_alleles' ",
                      "must have names on it"))
        what <- "some of the names on 'gapped_V_alleles'"
    } else {
        stop(wmsg("'gapped_V_alleles' must be a DNAStringSet object ",
                  "containing germline V gene allele gapped sequences, ",
                  "or the path to a FASTA file containing such sequences"))
    }
    names(gapped_V_alleles) <- clean_imgt_fasta_headers(allele_names, what)
    gapped_V_alleles
}

### Not exported!
### 'fwrcdr_widths' can be a named numeric vector or a numeric matrix with
### dimnames. Returns a named numeric vector or a numeric matrix with dimnames.
normarg_fwrcdr_widths <- function(fwrcdr_widths)
{
    if (!is.numeric(fwrcdr_widths))
        stop(wmsg("'fwrcdr_widths' must be an integer vector or matrix"))
    expected_len <- length(IMGT_DEFAULT_FWRCDR_WIDTHS)
    expected_nms <- names(IMGT_DEFAULT_FWRCDR_WIDTHS)
    expected_nms_in1string <- paste(expected_nms, collapse=", ")
    if (is.matrix(fwrcdr_widths)) {
        if (ncol(fwrcdr_widths) != expected_len)
            stop(wmsg("'fwrcdr_widths' must have ", expected_len, " columns"))
        dn <- dimnames(fwrcdr_widths)
        if (is.null(dn))
            stop(wmsg("'fwrcdr_widths' must have dimnames"))
        dn1 <- dn[[1L]]
        dn2 <- dn[[2L]]
        if (is.null(dn1) || is.null(dn2))
            stop(wmsg("'fwrcdr_widths' must have rownames and colnames"))
        if (!identical(dn2, expected_nms))
            stop(wmsg("the colnames on 'fwrcdr_widths' must be: ",
                      expected_nms_in1string))
        valid_rownames <- c(IG_LOCI, TR_LOCI)
        if (!all(dn1 %in% valid_rownames)) {
            in1string <- paste(valid_rownames, collapse=", ")
            stop(wmsg("valid rownames for 'fwrcdr_widths' are: ", in1string))
        }
        if (anyDuplicated(dn1))
            stop(wmsg("the rownames on 'fwrcdr_widths' cannot ",
                      "contain duplicates"))
    } else {
        if (length(fwrcdr_widths) != expected_len)
            stop(wmsg("'fwrcdr_widths' must have ", expected_len, " elements"))
        nms <- names(fwrcdr_widths)
        if (is.null(nms))
            stop(wmsg("'fwrcdr_widths' must have names"))
        if (!identical(nms, expected_nms))
            stop(wmsg("the names on 'fwrcdr_widths' must be: ",
                      expected_nms_in1string))
    }
    if (!is.integer(fwrcdr_widths))
        storage.mode(fwrcdr_widths) <- "integer"
    if (anyNA(fwrcdr_widths))
        stop(wmsg("'fwrcdr_widths' cannot contain NAs"))
    if (!all(fwrcdr_widths >= 1L))
        stop(wmsg("all values in 'fwrcdr_widths' must be >= 1"))
    fwrcdr_widths
}

### Not exported!
warn_if_allele_sequences_have_no_gaps <- function(ngaps)
{
    stopifnot(is.integer(ngaps))
    bad_idx <- which(ngaps == 0L)
    if (length(bad_idx) == 0L)
        return()
    allele_names <- names(ngaps)
    stopifnot(!is.null(allele_names))
    first_bad_allele <- allele_names[[bad_idx[[1L]]]]
    warning(wmsg(length(bad_idx), "/", length(ngaps), " V allele sequences ",
                 "have no gaps (e.g. allele ", first_bad_allele, ")"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_V_gene_delineations()
###

### Not exported!
EXTENDED_V_GENE_DELINEATION_COLNAMES <- c(
    "allele_name",
    V_GENE_DELINEATION_COLNAMES,
    "seq_len",
    "coding_frame_start",
    "starting_gap",
    "all_gaps_in_frame",
    "all_gaps_contained"
)

### Sanity checks.
stopifnot(
    identical(
        setdiff(names(NDM_DATA_COL2CLASS),
                EXTENDED_V_GENE_DELINEATION_COLNAMES),
        "chain_type"
    ),
    identical(
        setdiff(EXTENDED_V_GENE_DELINEATION_COLNAMES,
                names(NDM_DATA_COL2CLASS)),
        c("seq_len", "starting_gap", "all_gaps_in_frame", "all_gaps_contained")
    )
)

### Do the supplied ranges align with the underlying coding frame?
### Returns a logical vector parallel to 'dna_ranges'.
.dna_ranges_align_with_coding_frame <- function(dna_ranges)
{
    stopifnot(is(dna_ranges, "IRanges"))
    (end(dna_ranges) %% 3L == 0L) & (width(dna_ranges) %% 3L == 0L)
}

### Not used at the moment.
### All the ranges in 'dna_ranges' must align with the underlying coding
### frame that starts at position 1. An error will be raised if they don't.
.from_dna_to_aa_ranges <- function(dna_ranges)
{
    nms <- names(dna_ranges)
    if (is(dna_ranges, "PartitioningByWidth")) {
        widths <- setNames(width(dna_ranges), nms)
        ## Check alignment with underlying coding frame.
        stopifnot(all(widths %% 3L == 0L))
        return(PartitioningByWidth(widths %/% 3L))
    }
    if (is(dna_ranges, "PartitioningByEnd")) {
        ends <- setNames(end(dna_ranges), nms)
        ## Check alignment with underlying coding frame.
        stopifnot(all(ends %% 3L == 0L))
        return(PartitioningByEnd(ends %/% 3L))
    }
    if (!is(dna_ranges, "IntegerRanges"))
        stop(wmsg("'dna_ranges' must be an IRanges object ",
                  "or other IntegerRanges derivative"))
    stopifnot(all(.dna_ranges_align_with_coding_frame(dna_ranges)))
    ans_end <- end(dna_ranges) %/% 3L
    ans_width <- width(dna_ranges) %/% 3L
    IRanges(end=ans_end, width=ans_width, names=nms)
}

### Returns an integer vector parallel to 'imgt_bins' with the following
### attributes:
### - starting_gap: size (in number of nucleotides) of gap block located
###   at the very beginning of the sequence if any (0 if no such block);
### - all_gaps_in_frame: indicates whether the gap blocks align with the
###   underlying coding frame or not;
### - all_gaps_contained: TRUE if the gap blocks don't cross the FWR/CDR
###   boundaries.
.compute_fwrcdr_real_widths <- function(gap_pos, imgt_bins)
{
    stopifnot(is(gap_pos, "IRanges"), is(imgt_bins, "PartitioningByWidth"))
    gap_pos <- as(gap_pos, "StitchedIPos")

    ## Fastest way to check that the gap positions are strictly sorted.
    gap_blocks <- gap_pos@pos_runs  # IRanges object
    stopifnot(isNormal(gap_blocks))

    ## Note that some gapped V allele sequences from IMGT (e.g. human
    ## IGHV1-69-2*02 and IGHV7-34-1*01) can start with a gap block, and most
    ## of the time this block does not align with the underlying coding frame.
    if (length(gap_blocks) >= 1L && start(gap_blocks)[[1L]] == 1L) {
        starting_gap <- width(gap_blocks)[[1L]]
    } else {
        starting_gap <- 0L
    }

    ## Do all gap blocks align with the underlying coding frame?
    all_gaps_in_frame <- all(.dna_ranges_align_with_coding_frame(gap_blocks))

    ## Check that no gap block crosses FWR/CDR boundaries.
    counts <- countOverlaps(gap_blocks, imgt_bins, type="within")
    all_gaps_contained <- all(counts == 1L)

    imgt_bin_widths <- setNames(width(imgt_bins), names(imgt_bins))
    ngap_per_bin <- count_bin_hits(pos(gap_pos), imgt_bin_widths)
    ans <- imgt_bin_widths - ngap_per_bin

    attr(ans, "starting_gap") <- starting_gap
    attr(ans, "all_gaps_in_frame") <- all_gaps_in_frame
    attr(ans, "all_gaps_contained") <- all_gaps_contained
    ans
}

### Returns an ordinary data.frame.
.IRL_to_data_frame <- function(IRL)
{
    stopifnot(is(IRL, "CompressedIRangesList"))
    IRL_len <- length(IRL)
    all_ranges <- unlist(IRL, use.names=FALSE)
    expected_names <- rep.int(names(IMGT_DEFAULT_FWRCDR_WIDTHS), IRL_len)
    stopifnot(identical(expected_names, names(all_ranges)))

    idx0 <- seq_len(IRL_len) * length(IMGT_DEFAULT_FWRCDR_WIDTHS)
    df <- data.frame(
        allele_name=names(IRL),
        fwr1_start =start(all_ranges)[idx0 - 4L],
        fwr1_end   =end  (all_ranges)[idx0 - 4L],
        cdr1_start =start(all_ranges)[idx0 - 3L],
        cdr1_end   =end  (all_ranges)[idx0 - 3L],
        fwr2_start =start(all_ranges)[idx0 - 2L],
        fwr2_end   =end  (all_ranges)[idx0 - 2L],
        cdr2_start =start(all_ranges)[idx0 - 1L],
        cdr2_end   =end  (all_ranges)[idx0 - 1L],
        fwr3_start =start(all_ranges)[idx0],
        fwr3_end   =end  (all_ranges)[idx0]
    )
    df <- cbind(df, mcols(IRL, use.names=FALSE))  # ordinary data.frame
    stopifnot(identical(colnames(df), EXTENDED_V_GENE_DELINEATION_COLNAMES))
    df
}

.do_compute_V_gene_delineations <-
    function(gapped_V_alleles, fwrcdr_widths=IMGT_DEFAULT_FWRCDR_WIDTHS)
{
    stopifnot(is.integer(fwrcdr_widths),
              identical(names(fwrcdr_widths),
                        names(IMGT_DEFAULT_FWRCDR_WIDTHS)))
    ## IMGT FWR/CDR fixed intervals in nucleotide space.
    imgt_bins <- PartitioningByWidth(fwrcdr_widths * 3L)
    midx <- vmatchPattern(GAP_LETTER, gapped_V_alleles)

    ## Note that lengths() should propagate the names by default but it
    ## fails to do so on ByPos_MIndex object 'midx' at the moment (Biostrings
    ## 2.79.4).
    ## TODO: Fix this in Biostrings. Simplest fix is to define a lengths()
    ## method for MIndex objects that does 'lengths(endIndex(midx))' and
    ## add the names to that if 'use.names' is TRUE.
    ngaps <- setNames(lengths(midx), names(midx))
    warn_if_allele_sequences_have_no_gaps(ngaps)
    seq_len <- width(gapped_V_alleles) - ngaps  # lengths of ungapped sequences

    all_real_widths <- lapply(midx, .compute_fwrcdr_real_widths, imgt_bins)
    tmp <- lapply(all_real_widths, PartitioningByWidth)
    IRL <- as(tmp, "CompressedIRangesList")
    starting_gap <-
        vapply(all_real_widths, attr, integer(1), "starting_gap")
    all_gaps_in_frame <-
        vapply(all_real_widths, attr, logical(1), "all_gaps_in_frame")
    all_gaps_contained <-
        vapply(all_real_widths, attr, logical(1), "all_gaps_contained")
    coding_frame_start <- 2L - (starting_gap + 2L) %% 3L
    mcols(IRL) <- DataFrame(seq_len=seq_len,
                            coding_frame_start=coding_frame_start,
                            starting_gap=starting_gap,
                            all_gaps_in_frame=all_gaps_in_frame,
                            all_gaps_contained=all_gaps_contained)
    IRL
}

### 'gapped_V_alleles' can be a named DNAStringSet or BStringSet object,
### or the path to a FASTA file. Note that **all** the sequences
### in 'gapped_V_alleles' are expected to have gaps (the function will
### issue a warning if that's not the case).
### Returns a data.frame with 1 row per sequence in 'gapped_V_alleles'.
compute_V_gene_delineations <-
    function(gapped_V_alleles, fwrcdr_widths=IMGT_DEFAULT_FWRCDR_WIDTHS,
             as.IRangesList=FALSE)
{
    gapped_V_alleles <- normarg_gapped_V_alleles(gapped_V_alleles)
    fwrcdr_widths <- normarg_fwrcdr_widths(fwrcdr_widths)
    if (!isTRUEorFALSE(as.IRangesList))
        stop(wmsg("'as.IRangesList' must be TRUE or FALSE"))
    if (is.matrix(fwrcdr_widths)) {
        allele_loci <- substr(names(gapped_V_alleles), 1L, 3L)
        alleles_by_loci <- split(gapped_V_alleles, allele_loci)
        IRLs <- lapply(seq_along(alleles_by_loci),
            function(i) {
                locus <- names(alleles_by_loci)[[i]]
                if (locus %in% rownames(fwrcdr_widths)) {
                    locus_fwrcdr_widths <- fwrcdr_widths[locus, ]
                } else {
                    locus_fwrcdr_widths <- IMGT_DEFAULT_FWRCDR_WIDTHS
                }
                .do_compute_V_gene_delineations(alleles_by_loci[[i]],
                                   fwrcdr_widths=locus_fwrcdr_widths)
        })
        IRL <- unsplit2(IRLs, allele_loci)
    } else {
        IRL <- .do_compute_V_gene_delineations(gapped_V_alleles,
                                               fwrcdr_widths=fwrcdr_widths)
    }
    if (as.IRangesList)
        return(IRL)
    .IRL_to_data_frame(IRL)
}

compute_imgt_intdata <- function(...)
{
    .Deprecated("compute_V_gene_delineations")
    compute_V_gene_delineations(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_discordant_intdata()
### complete_intdata_with_ref()
###

### Not exported!
### See find_discordant_rows() in R/utils.R for what it returns.
find_discordant_intdata <- function(intdata, ref_intdata)
{
    check_ndm_data_col2class(intdata)
    check_ndm_data_col2class(ref_intdata)
    find_discordant_rows(intdata, ref_intdata, "allele_name")
}

### Not exported!
### See complete_df_with_ref() in R/utils.R for what it returns.
complete_intdata_with_ref <- function(intdata, ref_intdata,
                                      intdata_label="computed",
                                      ref_label="reference")
{
    check_ndm_data_col2class(intdata)
    check_ndm_data_col2class(ref_intdata)
    complete_df_with_ref(intdata, ref_intdata, "allele_name",
                         "\"internal data\"", intdata_label, ref_label)
}


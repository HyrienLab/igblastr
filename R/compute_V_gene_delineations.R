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

### There are some exceptions to this rule though: some organisms/loci have
### additional insertions in their FWR/CDR regions. IMGT online "Protein
### displays" service at https://www.imgt.org/IMGTrepertoire/Proteins/ shows
### these additional insertions.

### Not exported!
IMGT_MOUSE_FWRCDR_WIDTHS <- rbind(
    ## 1 additional insertion in FWR1 and 1 in FWR3 (between pos 84 & 85). See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=house%20mouse&latin=Mus%20musculus&group=TRAV
    TRA=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L),
    ## 1 additional insertion in FWR1 and 1 in FWR3 (between pos 84 & 85). See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=house%20mouse&latin=Mus%20musculus&group=TRDV
    TRD=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L)
)

### Not exported!
IMGT_RAT_FWRCDR_WIDTHS <- rbind(
    ## IMGT "Protein display" service shows 2 additional insertions in CDR1
    ## and 1 in FWR3. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Norway%20rat&latin=Rattus%20norvegicus&group=IGHV
    ## However, this seems very outdated! I sent an email to the IMGT team with
    ## subject "Protein display for rat IGHV seems out-of-sync with IGHV.fasta
    ## for rat" on Sep 11, 2026 about this.
    ## What we actually see in recent IGHV.fasta for rat is 1 additional
    ## insertion in CDR1 and 1 in FWR3.
    IGH=c(fwr1=26L, cdr1=13L, fwr2=17L, cdr2=10L, fwr3=40L)
)

### Not exported!
IMGT_RHESUS_MONKEY_FWRCDR_WIDTHS <- rbind(
    ## 2 additional insertions in FWR1 (between pos 15 & 16, and between
    ## pos 26 & 27, this 2nd insertion is considered part of the FWR1). See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=IGHV
    IGH=c(fwr1=28L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L),
    ## 1 additional insertion in FWR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=IGKV
    IGK=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L),
    ## 1 additional insertion in FWR1 and 2 in FWR2. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=IGLV
    IGL=c(fwr1=27L, cdr1=12L, fwr2=19L, cdr2=10L, fwr3=39L),
    ## 1 additional insertion in FWR1 and 1 in FWR3. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=TRAV
    TRA=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L),
    ## 1 additional insertion in FWR3. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=TRBV
    TRB=c(fwr1=26L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L),
    ## 1 additional insertion in FWR2. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=TRGV
    TRG=c(fwr1=26L, cdr1=12L, fwr2=18L, cdr2=10L, fwr3=39L)
)

### Not exported!
IMGT_CRAB_EATING_MACAQUE_FWRCDR_WIDTHS <- rbind(
    ## 1 additional insertion in FWR1 (between pos 26 & 27, considered part
    ## of the FWR1). See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=crab-eating%20macaque&latin=Macaca%20fascicularis&group=IGHV
    IGH=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L)
)

### Not exported!
IMGT_GORILLA_FWRCDR_WIDTHS <- rbind(
    ## 1 additional insertion in FWR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=western%20lowland%20gorilla&latin=Gorilla%20gorilla%20gorilla&group=TRAV
    TRA=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L),
    ## 1 additional insertion in FWR2. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=western%20lowland%20gorilla&latin=Gorilla%20gorilla%20gorilla&group=TRGV
    TRG=c(fwr1=26L, cdr1=12L, fwr2=18L, cdr2=10L, fwr3=39L)
)

### Not exported!
IMGT_CAT_FWRCDR_WIDTHS <- rbind(
    ## 1 additional insertion in FWR2 (between pos 46 & 47). See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20cat&latin=Felis%20catus&group=TRAV
    TRA=c(fwr1=26L, cdr1=12L, fwr2=18L, cdr2=10L, fwr3=39L)
)

### Not exported!
IMGT_DOG_FWRCDR_WIDTHS <- rbind(
    ## 2 additional insertions in CDR1, 1 in FWR2, and 1 in FWR3. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=dog&latin=Canis%20lupus%20familiaris&group=IGLV
    IGL=c(fwr1=26L, cdr1=14L, fwr2=18L, cdr2=10L, fwr3=40L),
    ## 1 additional insertion in FWR2 (between pos 46 & 47). See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=dog&latin=Canis%20lupus%20familiaris&group=TRAV
    TRA=c(fwr1=26L, cdr1=12L, fwr2=18L, cdr2=10L, fwr3=39L),
    ## 2 additional insertions in FWR3. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=dog&latin=Canis%20lupus%20familiaris&group=TRBV
    TRB=c(fwr1=26L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=41L)
)

### Not exported!
IMGT_BORNEAN_ORANGUTAN_FWRCDR_WIDTHS <- rbind(
    ## 3 additional insertions in FWR1 and 1 in FWR3. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Bornean%20orangutan&latin=Pongo%20pygmaeus&group=IGHV
    IGH=c(fwr1=29L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=40L)
)

### Not exported!
IMGT_COW_FWRCDR_WIDTHS <- rbind(
    ## 1 additional insertion in FWR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=bovine&latin=Bos%20taurus&group=TRBV
    TRB=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L),
    ## 2 additional insertions in FWR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=bovine&latin=Bos%20taurus&group=TRGV
    TRG=c(fwr1=28L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L),
    ## 6 additional insertions in CDR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=bovine&latin=Bos%20taurus&group=TRDV
    TRD=c(fwr1=26L, cdr1=18L, fwr2=17L, cdr2=10L, fwr3=39L)
)

### Not exported!
IMGT_PLATYPUS_FWRCDR_WIDTHS <- rbind(
    ## 1 additional insertion in CDR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=platypus&latin=Ornithorhynchus%20anatinus&group=IGHV
    IGH=c(fwr1=26L, cdr1=13L, fwr2=17L, cdr2=10L, fwr3=39L)
)

### Not exported!
IMGT_FERRET_FWRCDR_WIDTHS <- rbind(
    ## 1 additional insertion in FWR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20ferret&latin=Mustela%20putorius%20furo&group=IGKV
    IGK=c(fwr1=27L, cdr1=12L, fwr2=17L, cdr2=10L, fwr3=39L),
    ## 2 additional insertions in CDR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20ferret&latin=Mustela%20putorius%20furo&group=IGLV
    IGL=c(fwr1=26L, cdr1=14L, fwr2=17L, cdr2=10L, fwr3=39L),
    ## 1 additional insertion in FWR2 (between pos 46 & 47). See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20ferret&latin=Mustela%20putorius%20furo&group=TRAV
    TRA=c(fwr1=26L, cdr1=12L, fwr2=18L, cdr2=10L, fwr3=39L)
)

### Not exported!
IMGT_AMERICAN_MINK_FWRCDR_WIDTHS <-
    ## 1 additional insertion in FWR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=American%20mink&latin=Neogale%20vison&group=IGHV
    c(fwr1=27L, IMGT_DEFAULT_FWRCDR_WIDTHS[-1L])

### Not exported!
IMGT_ATLANTIC_SALMON_FWRCDR_WIDTHS <-
    ## 2 additional insertions in FWR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Atlantic%20salmon&latin=Salmo%20salar&group=IGHV
    c(fwr1=28L, IMGT_DEFAULT_FWRCDR_WIDTHS[-1L])

### Not exported!
IMGT_RAINBOW_TROUT_FWRCDR_WIDTHS <-
    ## 3 additional insertions in FWR1. See
    ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=rainbow%20trout&latin=Oncorhynchus%20mykiss&group=IGHV
    c(fwr1=29L, IMGT_DEFAULT_FWRCDR_WIDTHS[-1L])


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


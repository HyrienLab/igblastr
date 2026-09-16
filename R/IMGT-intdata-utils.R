### =========================================================================
### Low-level utilities to help compute internal data for IMGT organisms
### -------------------------------------------------------------------------
###


### Additional insertions in IMGT IG[HKL]V.fasta or TR[ABGD]V.fasta
### files can be seen by going to the "Protein displays" section at
### https://www.imgt.org/IMGTrepertoire/Proteins/ and clicking on any
### organism for any IG[HKL]V or TR[ABGD]V group.
### The '.IMGT_IG_ADDITIONAL_INSERTIONS' and '.IMGT_TR_ADDITIONAL_INSERTIONS'
### named lists below summarize these additional insertions for most (but not
### all) IMGT organisms. The two lists summarize the additional insertions in
### the sequences of the IG[HKL]V and TR[ABGD]V groups, respectively.
### Each list element must be either:
###   (1) a NULL;
###   (2) a named integer vector where the names are a subset
###       of 'names(IMGT_DEFAULT_FWRCDR_WIDTHS)';
###   (3) a named list where the names are IMGT V group names and each
###       list element is either a NULL or a named integer vector like (2).

.IMGT_IG_ADDITIONAL_INSERTIONS <- list(

    Bos_taurus=NULL,  # no additional insertions

    Canis_lupus_familiaris=list(
        ## No additional insertions in the IG[HK]V alleles.
        IGHV=NULL, IGKV=NULL,
        ## 2 additional insertions in CDR1, 1 in FWR2, and 1 in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=dog&latin=Canis%20lupus%20familiaris&group=IGLV
        IGLV=c(cdr1=2L, fwr2=1L, fwr3=1L)
    ),

    Ctenopharyngodon_idella=list(
        ## 1 additional insertion in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=grass%20carp&latin=Ctenopharyngodon%20idella&group=IGHV
        IGHV=c(fwr3=1L)
    ),

    Equus_caballus=list(
        ## 1 additional insertion in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=horse&latin=Equus%20caballus&group=IGHV
        IGHV=c(fwr3=1L),
        IGKV=NULL   # no additional insertions
    ),

    Felis_catus=list(
        ## 1 additional insertion in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20cat&latin=Felis%20catus&group=IGHV
        IGHV=c(fwr1=1L),
        IGKV=NULL,  # no additional insertions
        IGLV=NULL   # no additional insertions
    ),

    Gallus_gallus=list(
        IGHV=NULL,  # no additional insertions
        ## 1 additional insertion in FWR1 and 1 in FWR2. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=chicken&latin=Gallus%20gallus&group=IGLV
        IGLV=c(fwr1=1L, fwr2=1L)
    ),

    Gorilla_gorilla_gorilla=NULL,  # no additional insertions

    Homo_sapiens=NULL,  # no additional insertions

    Lemur_catta=NULL,  # no additional insertions

    Macaca_fascicularis=list(
        ## 1 additional insertion in FWR1 (between pos 26 & 27, considered
        ## part of the FWR1). See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=crab-eating%20macaque&latin=Macaca%20fascicularis&group=IGHV
        IGHV=c(fwr1=1L)
    ),

    Macaca_mulatta=list(
        ## 2 additional insertions in FWR1 (between pos 15 & 16, and between
        ## pos 26 & 27, this 2nd insertion is considered part of the FWR1). See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=IGHV
        IGHV=c(fwr1=2L),
        ## 1 additional insertion in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=IGKV
        IGKV=c(fwr1=1L),
        ## 1 additional insertion in FWR1 and 2 in FWR2. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=IGLV
        IGLV=c(fwr1=1L, fwr2=2L)
    ),

    Mus_musculus=NULL,  # no additional insertions

    Mustela_putorius_furo=list(
        IGHV=NULL,  # no additional insertions
        ## 1 additional insertion in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20ferret&latin=Mustela%20putorius%20furo&group=IGKV
        IGKV=c(fwr1=1L),
        ## 2 additional insertions in CDR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20ferret&latin=Mustela%20putorius%20furo&group=IGLV
        IGLV=c(cdr1=2L)
    ),

    Neogale_vison=list(
        ## 1 additional insertion in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=American%20mink&latin=Neogale%20vison&group=IGHV
        IGHV=c(fwr1=1L)
    ),

    Oncorhynchus_mykiss=list(
        ## 3 additional insertions in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=rainbow%20trout&latin=Oncorhynchus%20mykiss&group=IGHV
        IGHV=c(fwr1=3L)
    ),

    Ornithorhynchus_anatinus=list(
        ## 1 additional insertion in CDR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=platypus&latin=Ornithorhynchus%20anatinus&group=IGHV
        IGHV=c(cdr1=1L)
    ),

    Oryctolagus_cuniculus=NULL,  # no additional insertions

    Pongo_pygmaeus=list(
        ## 3 additional insertions in FWR1 and 1 in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Bornean%20orangutan&latin=Pongo%20pygmaeus&group=IGHV
        IGHV=c(fwr1=3L, fwr3=1L),
        IGKV=NULL,  # no additional insertions
        IGLV=NULL   # no additional insertions
    ),

    Rattus_norvegicus=list(
        ## IMGT "Protein display" service shows 2 additional insertions in CDR1
        ## and 1 in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Norway%20rat&latin=Rattus%20norvegicus&group=IGHV
        ## However, this seems very outdated! I sent an email to the IMGT team
        ## with subject "Protein display for rat IGHV seems out-of-sync with
        ## IGHV.fasta for rat" on Sep 11, 2026 about this.
        ## What we actually see in recent IGHV.fasta for rat is 1 additional
        ## insertion in CDR1 and 1 in FWR3.
        IGHV=c(cdr1=1L, fwr3=1L),
        IGKV=NULL,  # no additional insertions
        IGLV=NULL   # no additional insertions
    ),

    Salmo_salar=list(
        ## 2 additional insertions in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Atlantic%20salmon&latin=Salmo%20salar&group=IGHV
        IGHV=c(fwr1=2L)
    ),

    Sus_scrofa=NULL,  # no additional insertions

    Vicugna_pacos=list(
        IGHV=NULL   # no additional insertions
    )
)

.IMGT_TR_ADDITIONAL_INSERTIONS <- list(

    Bos_taurus=list(
        ## No additional insertions in the TRAV alleles.
        TRAV=NULL,
        ## 1 additional insertion in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=bovine&latin=Bos%20taurus&group=TRBV
        TRBV=c(fwr1=1L),
        ## 2 additional insertions in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=bovine&latin=Bos%20taurus&group=TRGV
        TRGV=c(fwr1=2L),
        ## 6 additional insertions in CDR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=bovine&latin=Bos%20taurus&group=TRDV
        TRDV=c(cdr1=6L)
    ),

    Camelus_dromedarius=NULL,  # no additional insertions

    Canis_lupus_familiaris=list(
        ## 1 additional insertion in FWR2 (between pos 46 & 47). See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=dog&latin=Canis%20lupus%20familiaris&group=TRAV
        TRAV=c(fwr2=1L),
        ## 2 additional insertions in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=dog&latin=Canis%20lupus%20familiaris&group=TRBV
        TRBV=c(fwr3=2L),
        TRGV=NULL,  # no additional insertions
        TRDV=NULL   # no additional insertions
    ),

    Danio_rerio=list(
        ## 1 additional insertion in FWR1, 1 in FWR2 (between pos 46 & 47),
        ## 2 in CDR2, and 2 in FWR3 (between pos 74 & 75 and between
        ## pos 84 & 85). See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=zebrafish&latin=Danio%20rerio&group=TRAV
        TRAV=c(fwr1=1L, fwr2=1L, cdr2=2L, fwr3=2L)
    ),

    Felis_catus=list(
        ## 1 additional insertion in FWR2 (between pos 46 & 47). See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20cat&latin=Felis%20catus&group=TRAV
        TRAV=c(fwr2=1L),
        ## No additional insertions in the TR[BGD]V alleles.
        TRBV=NULL, TRGV=NULL, TRDV=NULL
    ),

    Gorilla_gorilla_gorilla=list(
        ## 1 additional insertion in FWR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=western%20lowland%20gorilla&latin=Gorilla%20gorilla%20gorilla&group=TRAV
        TRAV=c(fwr1=1L),
        TRBV=NULL,  # no additional insertions
        ## 1 additional insertion in FWR2. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=western%20lowland%20gorilla&latin=Gorilla%20gorilla%20gorilla&group=TRGV
        TRGV=c(fwr2=1L),
        TRDV=NULL   # no additional insertions
    ),

    Heterocephalus_glaber=list(
        ## No additional insertions in the TR[ABG]V alleles.
        TRAV=NULL, TRBV=NULL, TRGV=NULL,
        ## 2 additional insertions in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=naked%20mole-rat&latin=Heterocephalus%20glaber&group=TRDV
        TRDV=c(fwr3=2L)
    ),

    Homo_sapiens=NULL,  # no additional insertions

    Macaca_fascicularis=list(
        TRBV=NULL   # no additional insertions
    ),

    Macaca_mulatta=list(
        ## 1 additional insertion in FWR1 and 1 in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=TRAV
        TRAV=c(fwr1=1L, fwr3=1L),
        ## 1 additional insertion in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=TRBV
        TRBV=c(fwr3=1L),
        ## 1 additional insertion in FWR2. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Rhesus%20monkey&latin=Macaca%20mulatta&group=TRGV
        TRGV=c(fwr2=1L),
        TRDV=NULL   # no additional insertions
    ),

    Mus_musculus=list(
        ## 1 additional insertion in FWR1 and 1 in FWR3 (between pos 84 & 85).
        ## See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=house%20mouse&latin=Mus%20musculus&group=TRAV
        TRAV=c(fwr1=1L, fwr3=1L),
        TRBV=NULL,  # no additional insertions
        TRGV=NULL,  # no additional insertions
        ## 1 additional insertion in FWR1 and 1 in FWR3 (between pos 84 & 85).
        ## See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=house%20mouse&latin=Mus%20musculus&group=TRDV
        TRDV=c(fwr1=1L, fwr3=1L)
    ),

    Mustela_putorius_furo=list(
        ## 1 additional insertion in FWR2 (between pos 46 & 47). See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=domestic%20ferret&latin=Mustela%20putorius%20furo&group=TRAV
        TRAV=c(fwr2=1L),
        ## No additional insertions in the TR[BGD]V alleles.
        TRBV=NULL, TRGV=NULL, TRDV=NULL
    ),

    Oryctolagus_cuniculus=NULL,  # no additional insertions

    Ovis_aries=list(
        ## 3 additional insertions in FWR2. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=sheep&latin=Ovis%20aries&group=TRAV
        TRAV=c(fwr2=3L),
        TRBV=NULL,  # no additional insertions
        ## 1 additional insertion in CDR1. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=sheep&latin=Ovis%20aries&group=TRDV
        TRDV=c(cdr1=1L)
    ),

    Pan_troglodytes=list(
        TRAV=NULL,  # no additional insertions
        TRBV=NULL,  # no additional insertions
        ## 1 additional insertion in FWR2. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=chimpanzee&latin=Pan%20troglodytes&group=TRGV
        TRGV=c(fwr2=1L),
        TRDV=NULL   # no additional insertions
    ),

    Pongo_abelii=NULL,  # no additional insertions

    Pongo_pygmaeus=list(
        TRBV=NULL,  # no additional insertions
        ## 1 additional insertion in FWR2. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=Bornean%20orangutan&latin=Pongo%20pygmaeus&group=TRGV
        TRGV=c(fwr2=1L)
    ),

    Sus_scrofa=list(
        ## 1 additional insertion in FWR2 (between pos 46 & 47). See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=pig&latin=Sus%20scrofa&group=TRBV
        TRBV=c(fwr2=1L),
        ## 1 additional insertion in FWR3. See
        ## https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=pig&latin=Sus%20scrofa&group=TRGV
        TRGV=c(fwr3=1L)
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_fwrcdr_widths_for_imgt_organism()
###

.make_vector_from_additional_insertions <- function(addins)
{
    ans <- IMGT_DEFAULT_FWRCDR_WIDTHS
    if (is.null(addins))
        return(ans)
    stopifnot(is.integer(addins))
    idx <- names(addins)
    stopifnot(!is.null(idx), all(idx %in% names(ans)))
    ans[idx] <- ans[idx] + addins
    ans
}

### Equivalent to:
###   do.call(rbind, lapply(addins, .make_vector_from_additional_insertions))
### but avoid the loop.
.make_matrix_from_additional_insertions <- function(addins, loci_prefix)
{
    stopifnot(is.list(addins))
    groups <- names(addins)
    valid_loci <- if (loci_prefix == "IG") IG_LOCI else TR_LOCI
    valid_groups <- paste0(valid_loci, "V")
    stopifnot(!is.null(groups), !anyDuplicated(groups),
              all(groups %in% valid_groups))
    ans <- matrix(IMGT_DEFAULT_FWRCDR_WIDTHS,
                  nrow=length(addins),
                  ncol=length(IMGT_DEFAULT_FWRCDR_WIDTHS), byrow=TRUE,
                  dimnames=list(groups, names(IMGT_DEFAULT_FWRCDR_WIDTHS)))
    unlisted_addins <- unlist(unname(addins), recursive=FALSE)
    if (is.null(unlisted_addins)) {
        unlisted_addins <- integer(0)
    } else {
        stopifnot(is.integer(unlisted_addins), !is.null(names(unlisted_addins)))
    }
    Midx <- cbind(rep.int(groups, lengths(addins)),
                  names(unlisted_addins))
    ans[Midx] <- ans[Midx] + unlisted_addins
    ans
}

get_fwrcdr_widths_for_imgt_organism <- function(organism,
                                                loci_prefix=c("IG", "TR"))
{
    organism <- normalize_IMGT_organism(organism)
    loci_prefix <- match.arg(loci_prefix)
    all_addins <- switch(loci_prefix,
                         IG=.IMGT_IG_ADDITIONAL_INSERTIONS,
                         TR=.IMGT_TR_ADDITIONAL_INSERTIONS,
                         stop(wmsg("invalid 'loci_prefix'")))
    if (!(organism %in% names(all_addins)))
        return(NA)
    addins <- all_addins[[organism]]
    if (is.list(addins))
        return(.make_matrix_from_additional_insertions(addins, loci_prefix))
    .make_vector_from_additional_insertions(addins)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### count_conserved_V_codons()
###

.get_amino_acids_at <- function(V_alleles, codon_ends)
{
    stopifnot(is(V_alleles, "DNAStringSet"), is.integer(codon_ends),
              length(V_alleles) == length(codon_ends))
    codon_starts <- codon_ends - 2L
    idx <- which(codon_starts >= 1L & codon_ends <= width(V_alleles))
    codons <- subseq(V_alleles[idx], codon_starts[idx], codon_ends[idx])
    translate(codons, no.init.codon=TRUE, if.fuzzy.codon="solve")
}

### Not exported!
### If the 2nd argument ('intdata') is missing then the 1st argument
### ('V_alleles') must be the name of an existing germline db.
### IMGT codons 23 (fwr1), 41 (fwr2), and 104 (fwr3) are expected to
### translate to C, W, and C, respectively (conserved codons).
### Returns a data.frame with 1 row per codon of interest and the
### following columns: C_count, W_count, other, percent_C, percent_W.
count_conserved_V_codons <- function(V_alleles, intdata)
{
    if (missing(intdata)) {
        check_germline_db_name(V_alleles)
        intdata <- load_intdata(V_alleles)
        V_alleles <- load_germline_sequences(V_alleles, region_types="V")
    }
    fwr1_ends   <- query_intdata(intdata, V_alleles, "fwr1_end", no.NAs=TRUE)
    fwr2_starts <- query_intdata(intdata, V_alleles, "fwr2_start", no.NAs=TRUE)
    fwr3_ends   <- query_intdata(intdata, V_alleles, "fwr3_end", no.NAs=TRUE)
    codons <- list(
        ## codon 23 + a few flanking codons:
        codon22 =.get_amino_acids_at(V_alleles, fwr1_ends - 12L),
        codon23 =.get_amino_acids_at(V_alleles, fwr1_ends - 9L),
        codon24 =.get_amino_acids_at(V_alleles, fwr1_ends - 6L),
        codon25 =.get_amino_acids_at(V_alleles, fwr1_ends - 3L),
        codon26 =.get_amino_acids_at(V_alleles, fwr1_ends),
        ## codon 41 + a few flanking codons:
        codon39 =.get_amino_acids_at(V_alleles, fwr2_starts + 2L),
        codon40 =.get_amino_acids_at(V_alleles, fwr2_starts + 5L),
        codon41 =.get_amino_acids_at(V_alleles, fwr2_starts + 8L),
        codon42 =.get_amino_acids_at(V_alleles, fwr2_starts + 11L),
        codon43 =.get_amino_acids_at(V_alleles, fwr2_starts + 14L),
        ## codon 104 + a few flanking codons:
        codon103=.get_amino_acids_at(V_alleles, fwr3_ends - 3L),
        codon104=.get_amino_acids_at(V_alleles, fwr3_ends),
        codon105=.get_amino_acids_at(V_alleles, fwr3_ends + 3L)
    )
    C_count <- vapply(codons, function(codon) sum(codon == "C"), integer(1))
    W_count <- vapply(codons, function(codon) sum(codon == "W"), integer(1))
    other <- vapply(codons, function(codon) sum(codon != "C" & codon != "W"),
                    integer(1))
    percent_C <- round(100 * C_count / lengths(codons), digits=2L)
    percent_W <- round(100 * W_count / lengths(codons), digits=2L)
    ans <- data.frame(C_count=C_count, W_count=W_count, other=other,
                      percent_C=percent_C, percent_W=percent_W)
    rownames(ans) <- names(codons)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_V_anchor_stats()
### compute_V_anchor_stats_for_germline_dbs_with_intdata()
###

### Not exported!
### If the 2nd argument ('intdata') is missing then the 1st argument
### ('V_alleles') must be the name of an existing germline db.
### Returns C23/W41/C104 percentages in a named numeric vector of length 3
### if 'add.stats.per.locus' is FALSE and in a 3-column matrix
### if 'add.stats.per.locus' is TRUE.
compute_V_anchor_stats <- function(V_alleles, intdata,
                                   add.stats.per.locus=FALSE)
{
    stopifnot(isTRUEorFALSE(add.stats.per.locus))
    if (missing(intdata)) {
        check_germline_db_name(V_alleles)
        intdata <- load_intdata(V_alleles)
        V_alleles <- load_germline_sequences(V_alleles, region_types="V")
    }
    conserved_codons <- c(codon23="C", codon41="W", codon104="C")
    Mindex <- cbind(names(conserved_codons),
                    paste0("percent_", conserved_codons))
    .compute_stats <- function(alleles) {
        stats <- count_conserved_V_codons(alleles, intdata)[Mindex]
        names(stats) <- paste0(conserved_codons,
                               sub("^codon", "", names(conserved_codons)))
        stats
    }
    if (!add.stats.per.locus)
        return(.compute_stats(V_alleles))
    allele_loci <- substr(names(V_alleles), 1L, 3L)
    ## The factor levels must be the unique allele loci in canonical order.
    f <- factor(allele_loci, levels=sort_unique_loci(allele_loci))
    alleles_by_locus <- split(V_alleles, f)
    do.call(rbind, c(list(all_loci=.compute_stats(V_alleles)),
                     lapply(alleles_by_locus, .compute_stats)))
}

### Returns a vector of db names parallel to 'loci'.
.form_single_locus_germline_db_names <- function(db_name, loci)
{
    stopifnot(isSingleNonWhiteString(db_name), is.character(loci))
    parts <- strsplit(db_name, ".", fixed=TRUE)[[1L]]
    loci_idx <- which(grepl("\\<IG[HKL]|TR[ABGD]\\>", parts))
    if (length(loci_idx) == 0L)
        return(paste0(db_name, ".", loci, recycle0=TRUE))
    n1 <- loci_idx[[1L]] - 1L
    prefix <- paste(paste0(head(parts, n=n1), ".", recycle0=TRUE), collapse="")
    n2 <- length(parts) - loci_idx[[1L]]
    suffix <- paste(paste0(".", tail(parts, n=n2), recycle0=TRUE), collapse="")
    paste0(prefix, loci, suffix, recycle0=TRUE)
}

### Not exported!
compute_V_anchor_stats_for_germline_dbs_with_intdata <-
    function(add.stats.per.locus=FALSE)
{
    all_db_names <- list_germline_dbs(with.intdata.only=TRUE, names.only=TRUE)
    all_stats <- lapply(setNames(all_db_names, all_db_names),
        function(db_name) {
            stats <- compute_V_anchor_stats(db_name,
                             add.stats.per.locus=add.stats.per.locus)
            if (!add.stats.per.locus) {
                stopifnot(is.vector(stats))
                return(stats)  # returns a named vector
            }
            ## Will return a matrix.
            stopifnot(is.matrix(stats), nrow(stats) >= 2L)
            rnms <- rownames(stats)
            stopifnot(identical(rnms[[1L]], "all_loci"))
            single_locus_db_names <-
                .form_single_locus_germline_db_names(db_name, rnms[-1L])
            rownames(stats) <- c(db_name, single_locus_db_names)
            if (nrow(stats) > 2L)
                return(stats)
            row1 <- stats[1L, , drop=FALSE]
            if (identical(row1, stats[2L, , drop=FALSE])) row1 else stats
        })
    do.call(rbind, all_stats)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Percent cysteines (C) at codon23 & codon104 and percent tryptophans (W)
### at codon41.

### Results obtained with igblastr 1.3.22 for IMGT 202631-1.

### Output of:
###   reset_germline_dbs()
###   igblastr:::install_all_IMGT_germline_dbs("202631-1")
###   igblastr:::compute_V_anchor_stats_for_germline_dbs_with_intdata()
###                                                          C23    W41   C104
### IMGT-202631-1.Bos_taurus.IGH+IGK+IGL                   98.88  95.51  95.00
### IMGT-202631-1.Canis_lupus_familiaris.IGH+IGK+IGL       97.89  95.36  96.20
### IMGT-202631-1.Ctenopharyngodon_idella.IGH              99.29  89.29 100.00
### IMGT-202631-1.Equus_caballus.IGH+IGK                   91.01  88.76  93.26
### IMGT-202631-1.Felis_catus.IGH+IGK+IGL                  91.47  96.90  93.02
### IMGT-202631-1.Gallus_gallus.IGH+IGL                    98.37  97.01  98.44
### IMGT-202631-1.Gorilla_gorilla_gorilla.IGH+IGK+IGL      98.26  89.58  93.40
### IMGT-202631-1.Homo_sapiens.IGH+IGK+IGL                 98.47  94.54  97.42
### IMGT-202631-1.Lemur_catta.IGH+IGK+IGL                  96.92  99.56 100.00
### IMGT-202631-1.Macaca_fascicularis.IGH                  93.18  95.45  98.85
### IMGT-202631-1.Macaca_mulatta.IGH+IGK+IGL               95.62  93.87  96.05
### IMGT-202631-1.Mus_musculus.IGH+IGK+IGL                 98.49  94.68  94.51
### IMGT-202631-1.Mustela_putorius_furo.IGH+IGK+IGL        95.48  93.22  94.35
### IMGT-202631-1.Neogale_vison.IGH                        96.77 100.00  96.77
### IMGT-202631-1.Oncorhynchus_mykiss.IGH                  98.73  92.36  94.90
### IMGT-202631-1.Ornithorhynchus_anatinus.IGH             97.78 100.00  95.45
### IMGT-202631-1.Oryctolagus_cuniculus.IGH+IGK+IGL        98.64 100.00 100.00
### IMGT-202631-1.Pongo_pygmaeus.IGH+IGK+IGL               95.17  94.18  93.47
### IMGT-202631-1.Rattus_norvegicus.IGH+IGK+IGL            95.76  91.56  95.48
### IMGT-202631-1.Salmo_salar.IGH                          95.33  92.00  91.33
### IMGT-202631-1.Sus_scrofa.IGH+IGK+IGL                   96.92  96.92  95.38
### IMGT-202631-1.Vicugna_pacos.IGH                        98.81  97.62 100.00

### Output of:
###   reset_germline_dbs()
###   igblastr:::install_all_IMGT_germline_dbs("202631-1", tcr.db=TRUE)
###   igblastr:::compute_V_anchor_stats_for_germline_dbs_with_intdata()
###                                                          C23    W41   C104
### IMGT-202631-1.Bos_taurus.TRA+TRB+TRG+TRD               98.16  94.23  99.21
### IMGT-202631-1.Camelus_dromedarius.TRA+TRB+TRG+TRD     100.00 100.00 100.00
### IMGT-202631-1.Canis_lupus_familiaris.TRA+TRB+TRG+TRD   98.77  98.77  98.77
### IMGT-202631-1.Danio_rerio.TRA+TRD                      99.64  98.93 100.00
### IMGT-202631-1.Felis_catus.TRA+TRB+TRG+TRD              98.85  98.85 100.00
### IMGT-202631-1.Gorilla_gorilla_gorilla.TRA+TRB+TRG+TRD  97.34  99.47 100.00
### IMGT-202631-1.Heterocephalus_glaber.TRA+TRB+TRG+TRD    95.65  95.65  96.52
### IMGT-202631-1.Homo_sapiens.TRA+TRB+TRG+TRD             96.94  98.90  97.75
### IMGT-202631-1.Macaca_fascicularis.TRB                  95.45  98.48  98.46
### IMGT-202631-1.Macaca_mulatta.TRA+TRB+TRG+TRD           97.69  99.54  98.60
### IMGT-202631-1.Mus_musculus.TRA+TRB+TRG+TRD             99.47  97.16  98.70
### IMGT-202631-1.Mustela_putorius_furo.TRA+TRB+TRG+TRD    96.64  94.12  97.48
### IMGT-202631-1.Oryctolagus_cuniculus.TRA+TRB+TRG+TRD    99.32  99.32 100.00
### IMGT-202631-1.Ovis_aries.TRA+TRB+TRD                   99.47  94.93  99.20
### IMGT-202631-1.Pan_troglodytes.TRA+TRB+TRG+TRD          97.59  98.19  97.59
### IMGT-202631-1.Pongo_abelii.TRA+TRB+TRG+TRD             96.15  96.92  85.77
### IMGT-202631-1.Pongo_pygmaeus.TRB+TRG                   95.73  96.58  96.58
### IMGT-202631-1.Sus_scrofa.TRB+TRG                       92.68  97.56 100.00


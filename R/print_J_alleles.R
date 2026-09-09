### =========================================================================
### print_J_alleles()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_Jparts()
###

.get_aa_parts <- function(J_aa, cdr3_coding_frame_width)
{
    cdr3_ncodon <- cdr3_coding_frame_width %/% 3L
    cdr3 <- subseq(J_aa, end=cdr3_ncodon)
    fwr4 <- subseq(J_aa, start=cdr3_ncodon+1L)
    DataFrame(cdr3=cdr3, fwr4=fwr4)
}

.get_dna_parts <- function(J_dna, cdr3_end)
{
    cdr3 <- subseq(J_dna, end=cdr3_end)
    fwr4 <- subseq(J_dna, start=cdr3_end+1L)
    DataFrame(cdr3=cdr3, fwr4=fwr4)
}

.extract_Jparts <- function(J_alleles, auxdata,
                            translate=FALSE, igblast_organism=NA)
{
    if (!is(J_alleles, "DNAStringSet"))
        stop(wmsg("'J_alleles' must be a DNAStringSet object"))
    fasta_headers <- names(J_alleles)
    if (is.null(fasta_headers))
        stop(wmsg("'J_alleles' must have names"))
    names(J_alleles) <- clean_imgt_fasta_headers(fasta_headers)
    if (!isTRUEorFALSE(translate))
        stop(wmsg("'translate' must be TRUE or FALSE"))

    coding_frame_start <-
            query_auxdata(auxdata, J_alleles, "coding_frame_start",
                          no.NAs=TRUE)
    cdr3_end <- query_auxdata(auxdata, J_alleles, "cdr3_end")  # 0-based
    splitidx <- which(!is.na(cdr3_end))
    cdr3_end <- cdr3_end + 1L  # 1-based
    cdr3_coding_frame_width <- cdr3_end - coding_frame_start
    stopifnot(all(cdr3_coding_frame_width %% 3L == 0L, na.rm=TRUE))
    if (translate) {
        full_seq <- translate_codons(J_alleles, offset=coding_frame_start)
        aa0 <- rep.int(AAStringSet(""), length(full_seq))
        parts <- DataFrame(cdr3=aa0, fwr4=aa0)
        parts[splitidx, ] <- .get_aa_parts(full_seq[splitidx],
                                           cdr3_coding_frame_width[splitidx])
    } else {
        full_seq <- J_alleles
        ## Drop 'full_seq' metadata cols so they don't do strange things
        ## when we pass 'full_seq' to the DataFrame() call below.
        mcols(full_seq) <- NULL
        dna0 <- rep.int(DNAStringSet(""), length(full_seq))
        parts <- DataFrame(cdr3=dna0, fwr4=dna0)
        parts[splitidx, ] <- .get_dna_parts(full_seq[splitidx],
                                            cdr3_end[splitidx])
    }
    allele_names <- names(J_alleles)
    ans <- cbind(DataFrame(allele_name=allele_names),
                 parts,
                 DataFrame(full_seq=full_seq))
    if (!identical(igblast_organism, NA)) {
        if (!isSingleNonWhiteString(igblast_organism))
            stop(wmsg("'igblast_organism' must be NA or ",
                      "a single (non-empty) string"))
        igblast_organism <- normalize_igblast_organism(igblast_organism)
        igblast_auxdata <- load_auxdata(igblast_organism)
        ans$also_in_IgBLAST_auxdata <-
                    allele_names %in% igblast_auxdata[ , "allele_name"]
    }
    extra_cols <- mcols(J_alleles, use.names=FALSE)
    if (length(extra_cols) != 0L)
        ans <- cbind(ans, extra_cols)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .print_Jparts()
###

.make_Jparts_header_line <- function(labelW, cdr3W, cdr3fwr4_sep, fwr4W,
                                     extra_colnames)
{
    stopifnot(isSingleInteger(labelW),
              isSingleInteger(cdr3W),
              isSingleString(cdr3fwr4_sep),
              isSingleInteger(fwr4W),
              is.null(extra_colnames) || is.character(extra_colnames))
    if (cdr3W >= 5L) {
        cdr3 <- paste0(strrep(" ", cdr3W-5L), "CDR3>")
    } else if (cdr3W == 4L) {
        cdr3 <- "CDR3"
    } else {
        cdr3 <- strrep(" ", cdr3W)
    }
    if (fwr4W >= 5L) {
        fwr4 <- paste0("<FWR4", strrep(" ", fwr4W-5L))
    } else if (fwr4W == 4L) {
        fwr4 <- "FWR4"
    } else {
        fwr4 <- strrep(" ", fwr4W)
    }
    ans <- paste0(strrep(" ", labelW), cdr3, cdr3fwr4_sep, fwr4)
    if (length(extra_colnames) != 0L)
        ans <- paste0(ans, paste(extra_colnames, collapse=" "))
    ans
}

.make_Jparts_line <- function(i, labels, Lfillers, cdr3,
                                 cdr3fwr4_sep, fwr4, Rfillers,
                                 full_seq, extra_cols)
{
    stopifnot(isSingleInteger(i),
              is.character(labels),
              is.character(Lfillers),
              is(cdr3, "XStringSet"),
              isSingleString(cdr3fwr4_sep),
              is(fwr4, "XStringSet"),
              is.character(Rfillers),
              is.matrix(extra_cols),
              length(labels) == length(Lfillers),
              length(labels) == length(cdr3),
              length(labels) == length(fwr4),
              length(labels) == length(Rfillers),
              length(labels) == nrow(extra_cols))
    seq1 <- paint_XString(cdr3[[i]])
    seq2 <- paint_XString(fwr4[[i]])
    fseq <- full_seq[[i]]
    ## Both 'Ltag' and 'Rtag' must be single characters. If that
    ## needs to change then the call to .make_Jparts_header_line()
    ## in .print_Jparts() below will need to be adjusted.
    if (length(fseq) == 0L) {
        Ltag <- Rtag <- " "
        stuff <- c(Lfillers[[i]], seq1, cdr3fwr4_sep, seq2, Rfillers[[i]])
    } else {
        Ltag <- Rtag <- "?"
        Linnertag <- ">"
        Rinnertag <- " <"  # note the space before the <
        LinnertagW <- nchar(Lfillers[[i]]) + nchar(cdr3fwr4_sep) +
                      nchar(Rfillers[[i]]) - nchar(fseq) - nchar(Rinnertag)
        Linnertag <- format(Linnertag, width=LinnertagW)
        stuff <- c(Linnertag, paint_XString(fseq), Rinnertag)
    }
    ans <- paste0(labels[[i]], Ltag, paste(stuff, collapse=""), Rtag)
    if (ncol(extra_cols) != 0L)
        ans <- paste0(ans, paste(extra_cols[i, ], collapse=" "))
    ans
}

.print_Jparts <- function(Jparts, filler=".", cdr3fwr4_sep=" ")
{
    stopifnot(is(Jparts, "DataFrame"))
    allele_names <- Jparts[ , "allele_name"]
    cdr3         <- Jparts[ , "cdr3"]
    fwr4         <- Jparts[ , "fwr4"]
    full_seq     <- Jparts[ , "full_seq"]
    splitidx <- which(width(cdr3) != 0L | width(fwr4) != 0L)
    full_seq[splitidx] <- ""
    core_colnames <- c("allele_name", "cdr3", "fwr4", "full_seq")
    core_col_idx <- match(core_colnames, colnames(Jparts))
    extra_cols <- Jparts[-core_col_idx]

    ## Format all columns.
    labels <- paste0(format(seq_len(nrow(Jparts))), ". ",
                     format(paste0(allele_names, ": ")))
    labelW <- nchar(labels[[1L]])
    stopifnot(all(nchar(labels) == labelW))
    cdr3_maxwidth <- max(width(cdr3))
    Lfillers <- strrep(filler, cdr3_maxwidth - width(cdr3))
    fwr4_maxwidth <- max(width(fwr4))
    Rfillers <- strrep(filler, fwr4_maxwidth - width(fwr4))
    extra_cols <- DataFrame_as_formatted_matrix(extra_cols)
    extra_colnames <- colnames(extra_cols)

    ## Print everything.

    locus <- substr(allele_names, 1L, 3L)
    group_lens <- runLength(Rle(locus))
    if (length(group_lens) == length(unique(locus))) {
        groupend_idx <- cumsum(head(group_lens, n=-1L))
    } else {
        groupend_idx <- integer(0)
    }

    ## We add 1 to 'cdr3_maxwidth' and 'fwr4_maxwidth' to account for
    ## the size of the 'Ltag' and 'Rtag' used in .make_Jparts_line().
    ## See .make_Jparts_line() above in this file.
    line <- .make_Jparts_header_line(labelW, cdr3_maxwidth+1L,
                                     cdr3fwr4_sep, fwr4_maxwidth+1L,
                                     extra_colnames)
    stopifnot(isSingleString(line))
    message(line)
    for (i in seq_len(nrow(Jparts))) {
        line <- .make_Jparts_line(i, labels, Lfillers, cdr3,
                                     cdr3fwr4_sep, fwr4, Rfillers,
                                     full_seq, extra_cols)
        stopifnot(isSingleString(line))
        message(line)
        if (i %in% groupend_idx)
            message("")
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### print_J_alleles()
###

### Displays the J alleles sequences justified with respect to their
### CDR3/FWR4 junction.
print_J_alleles <- function(J_alleles, auxdata, translate=FALSE,
                            igblast_organism=NA, filler=".", cdr3fwr4_sep=" ")
{
    if (!is(J_alleles, "DNAStringSet")) {
        if (!isSingleNonWhiteString(J_alleles))
            stop(wmsg("'J_alleles' must be a DNAStringSet object ",
                      "containing germline J gene allele sequences, ",
                      "or a single string that is the name of a cached ",
                      "germline db"))
        db_name <- J_alleles
        J_alleles <- load_germline_sequences(db_name, region_types="J")
        if (!missing(auxdata))
            stop(wmsg("'auxdata' should not be supplied when 'J_alleles' ",
                      "is the name of a cached germline db"))
        auxdata <- load_auxdata(db_name)
        if (identical(igblast_organism, NA)) {
            igblast_organism <- infer_igblast_organism_from_db_name(db_name)
            if (is.na(igblast_organism))
                igblast_organism <- NA  # replace NA_character_ with NA
        }
    }
    if (!isSingleString(filler) || nchar(filler) != 1L)
        stop(wmsg("'filler' must be a single character"))
    if (!isSingleString(cdr3fwr4_sep))
        stop(wmsg("'cdr3fwr4_sep' must be a single string"))
    Jparts <- .extract_Jparts(J_alleles, auxdata,
                              translate=translate,
                              igblast_organism=igblast_organism)
    .print_Jparts(Jparts, filler=filler, cdr3fwr4_sep=cdr3fwr4_sep)
    invisible(Jparts)
}


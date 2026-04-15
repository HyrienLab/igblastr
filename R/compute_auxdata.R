### =========================================================================
### Compute IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .find_heavy_fwr4_starts()
### .find_light_fwr4_starts()
###

### For alleles in the IGHJ group (i.e. BCR germline J gene alleles on the
### heavy chain), the FWR4 region is expected to start with AA motif "WGXG.
.WGXG_pattern <- "TGGGGNNNNGGN"  # reverse-translation of "WGXG"

### For all other J alleles, that is, for alleles in the IG[KL]J groups
### (i.e. BCR germline J gene alleles on the light chain) and all TCR
### germline J gene alleles, the FWR4 region is expected to start with
### AA motif "FGXG".
.FGXG_pattern <- "TTYGGNNNNGGN"  # reverse-translation of "FGXG"

### EXPERIMENTAL!
### The "FGXG" motif is not found for 4 J alleles in
### IMGT-202531-1.Mus_musculus.IGH+IGK+IGL: IGKJ3*01, IGKJ3*02,
### IGLJ2P*01, IGLJ3P*01. However, except for IGLJ2P*01, these alleles
### are annotated in mouse_gl.aux with a CDR3 end reported at position 6
### (0-based). Turns out that for the 3 alleles annotated in mouse_gl.aux,
### the two first codons of the FWR4 region translate to AA sequence "FS".
### Is this a coincidence or does the FS sequence actually play a role on
### the light chain? What do biologists say about this? In particular, does
### it make sense to use this alternative motif to identify the start of
### the FWR4 region on the light chain when the "FGXG" motif is not found?
### Note that all the possible reverse-translations of FS cannot be
### represented with a single DNA pattern (even with the use of IUPAC
### ambiguity codes).
.FS_pattern1 <- "TTYTCN"
.FS_pattern2 <- "TTYAGY"

### UPDATE on using the "FS" motif to identify the start of the FWR4
### region on the light chain when the "FGXG" motif is not found:
### Works well for IMGT-202531-1.Mus_musculus.IGH+IGK+IGL (well, it was
### specifically designed for that so no surprise here), but not
### so well for IMGT-202531-1.Rattus_norvegicus.IGH+IGK+IGL or
### IMGT-202531-1.Oryctolagus_cuniculus.IGH+IGK+IGL (rabbit)
### or IMGT-202531-1.Macaca_mulatta.IGH+IGK+IGL (rhesus monkey).
### So we disabled this feature in .find_light_fwr4_starts() below.

### .find_heavy_fwr4_starts() and .find_light_fwr4_starts() both return
### a named integer vector parallel to 'J_alleles' that contains
### the **0-based** FWR4 start position for each sequence in 'J_alleles'.
### Th FWR4 start will be set to NA for alleles that don't have a match.
### For alleles with more than one match, we keep the first match only.
### The names on the returned vector indicate the AA motif that was used
### to determine the start of the FWR4 region.

.find_heavy_fwr4_starts <- function(J_alleles)
{
    stopifnot(is(J_alleles, "DNAStringSet"))
    m <- vmatchPattern(.WGXG_pattern, J_alleles, fixed=FALSE)
    ans <- as.integer(heads(start(m), n=1L)) - 1L
    names(ans) <- ifelse(is.na(ans), NA_character_, "WGXG")
    ans
}

.find_light_fwr4_starts <- function(J_alleles)
{
    stopifnot(is(J_alleles, "DNAStringSet"))
    m <- vmatchPattern(.FGXG_pattern, J_alleles, fixed=FALSE)
    FGXG_starts <- as.integer(heads(start(m), n=1L))
    names(FGXG_starts) <- ifelse(is.na(FGXG_starts), NA_character_, "FGXG")
    ## Disabling search for alternative "FS" motif for now.
    #na_idx <- which(is.na(FGXG_starts))
    #if (length(na_idx) != 0L) {
    #    dangling_alleles <- J_alleles[na_idx]
    #    m <- vmatchPattern(.FS_pattern1, dangling_alleles, fixed=FALSE)
    #    FS_starts1 <- as.integer(heads(start(m), n=1L))
    #    m <- vmatchPattern(.FS_pattern2, dangling_alleles, fixed=FALSE)
    #    FS_starts2 <- as.integer(heads(start(m), n=1L))
    #    FS_starts <- pmin(FS_starts1, FS_starts2, na.rm=TRUE)
    #    FGXG_starts[na_idx] <- FS_starts
    #    names(FGXG_starts)[na_idx] <-
    #        ifelse(is.na(FS_starts), NA_character_, "FS")
    #}
    FGXG_starts - 1L
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_auxdata()
###

.VALID_J_GROUPS <- paste0("IG", c("H", "K", "L"), "J")

### Returns a data.frame with the same column names as the data.frame
### returned by load_auxdata() (see file R/auxdata-utils.R), minus
### the "extra_bps" column and plus the "fwr4_start_motif" column.
### NOTE: We set coding_frame_start/cdr3_end/fwr4_start_motif
### to NA for alleles for which the FWR4 start cannot be determined.
.compute_auxdata_for_J_group <- function(J_alleles, J_group)
{
    stopifnot(is(J_alleles, "DNAStringSet"), J_group %in% .VALID_J_GROUPS)
    allele_names <- names(J_alleles)
    stopifnot(!is.null(allele_names))
    if (length(J_alleles) == 0L) {
        chain_type <- character(0)
    } else {
        allele_groups <- substr(allele_names, 1L, 4L)
        stopifnot(all(allele_groups == J_group))
        chain_type <- paste0("J", substr(J_group, 3L, 3L))
    }
    if (J_group == "IGHJ") {
        fwr4_starts <- .find_heavy_fwr4_starts(J_alleles)
    } else {
        fwr4_starts <- .find_light_fwr4_starts(J_alleles)
    }
    coding_frame_starts <- fwr4_starts %% 3L
    data.frame(
        allele_name       =allele_names,
        coding_frame_start=coding_frame_starts,  # 0-based
        chain_type        =chain_type,
        cdr3_end          =fwr4_starts - 1L      # 0-based
        ## Returning this column only made sense when we were using "FS"
        ## motif as a 2nd-chance motif on the light chain.
        #fwr4_start_motif  =names(fwr4_starts)
    )
}

.normarg_codon_starts <- function(codon_starts, allele_names)
{
    if (!is.numeric(codon_starts))
        stop(wmsg("'codon_starts' must be NULL or an integer vector"))
    if (is.null(names(codon_starts)))
        stop(wmsg("'codon_starts' must carry the names of the J alleles"))
    if (!identical(names(codon_starts), allele_names))
        stop(wmsg("the names on 'J_alleles' and 'codon_starts' must ",
                  "be identical"))
    if (!is.integer(codon_starts))
        codon_starts <- setNames(as.integer(codon_starts), names(codon_starts))
    if (!(all(codon_starts %in% c(1:3, NA_integer_))))
        stop(wmsg("the non-NA values in 'codon_starts' must ",
                  "be >= 1 and <= 3"))
    codon_starts
}

.refine_with_codon_starts <- function(auxdata, codon_starts)
{
    stopifnot(is.data.frame(auxdata), is.integer(codon_starts),
              identical(auxdata[ , "allele_name"], names(codon_starts)))
    coding_frame_start <- auxdata[ , "coding_frame_start"]
    expected_coding_frame_start <- codon_starts - 1L
    bad_idx <- which(coding_frame_start != expected_coding_frame_start)
    if (length(bad_idx) != 0L) {
        in1string <- paste(auxdata[bad_idx, "allele_name"], collapse=", ")
        stop(wmsg("the supplied \"codon start\" is in disagreement with ",
                  "the computed \"coding frame start\" for allele(s): ",
                  in1string))
    }
    na_idx <- which(is.na(coding_frame_start))
    if (length(na_idx) != 0L) {
        coding_frame_start[na_idx] <- expected_coding_frame_start[na_idx]
        auxdata$coding_frame_start <- coding_frame_start
    }
    auxdata
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
        warning(wmsg("the \"coding frame start\" could not be determined ",
                     "for J allele(s): ", in1string),
                "\n  ",
                wmsg("--> coding_frame_start, cdr3_end, and extra_bps ",
                     "were set to NA for these alleles"))
    }
}

### Returns a data.frame with 1 row per sequence in 'J_alleles'.
compute_auxdata <- function(J_alleles, codon_starts=NULL, no.warning=FALSE)
{
    if (!is(J_alleles, "DNAStringSet"))
        stop(wmsg("'J_alleles' must be DNAStringSet object"))
    allele_names <- names(J_alleles)
    if (is.null(allele_names))
        stop(wmsg("'J_alleles' must have names"))
    names(J_alleles) <- allele_names <- clean_imgt_fasta_headers(allele_names)

    allele_groups <- substr(allele_names, 1L, 4L)
    if (!all(allele_groups %in% .VALID_J_GROUPS))
        stop(wmsg("all allele names must start with 'IG[HKL]J'"))

    if (!is.null(codon_starts))
        codon_starts <- .normarg_codon_starts(codon_starts, allele_names)

    if (!isTRUEorFALSE(no.warning))
        stop(wmsg("'no.warning' must be TRUE or FALSE"))

    JH_alleles <- J_alleles[allele_groups == "IGHJ"]
    JK_alleles <- J_alleles[allele_groups == "IGKJ"]
    JL_alleles <- J_alleles[allele_groups == "IGLJ"]
    JH_df <- .compute_auxdata_for_J_group(JH_alleles, "IGHJ")
    JK_df <- .compute_auxdata_for_J_group(JK_alleles, "IGKJ")
    JL_df <- .compute_auxdata_for_J_group(JL_alleles, "IGLJ")
    ans <- rbind(JH_df, JK_df, JL_df)
    rownames(ans) <- NULL

    i <- match(allele_names, ans[ , "allele_name"])
    ans <- S4Vectors:::extract_data_frame_rows(ans, i)

    if (!is.null(codon_starts))
        ans <- .refine_with_codon_starts(ans, codon_starts)

    ## Add "extra_bps" column.
    ans$extra_bps <- (width(J_alleles) - ans[ , "coding_frame_start"]) %% 3L

    if (!no.warning)
        .warn_if_missing_coding_frame_starts(ans)

    ans
}


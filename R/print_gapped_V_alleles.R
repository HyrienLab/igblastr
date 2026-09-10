### =========================================================================
### print_gapped_V_alleles()
### -------------------------------------------------------------------------
###


.VPART_NAMES <- paste0(c("fwr", "cdr"), rep(1:3, each=2))
stopifnot(
    identical(c(.VPART_NAMES, "fwr4"), FWRCDR_NAMES),
    identical(.VPART_NAMES, c(names(IMGT_DEFAULT_FWRCDR_WIDTHS), "cdr3"))
)

.normarg_region_types2 <- function(region_types=NULL)
{
    if (is.null(region_types))
        return(.VPART_NAMES)
    if (!is.character(region_types))
        stop(wmsg("'region_types' must be NULL or a character vector"))
    if (length(region_types) == 0L)
        stop(wmsg("'region_types' cannot be an empty character vector"))
    if (anyNA(region_types) || anyDuplicated(region_types))
        stop(wmsg("'region_types' cannot contain NAs or duplicates"))
    m <- match(region_types, .VPART_NAMES)
    if (anyNA(m) || is.unsorted(m, strictly=TRUE)) {
        in1string <- paste0("\"", .VPART_NAMES, "\"", collapse=", ")
        stop(wmsg("'region_types' must be an ordered subset ",
                  "of 'c(", in1string, ")'"))
    }
    region_types
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_Vparts()
###

.extract_Vparts <- function(gapped_V_alleles, fwrcdr_widths, region_types,
                            igblast_organism=NA)
{
    stopifnot(is(gapped_V_alleles, "XStringSet"),
              is.integer(fwrcdr_widths),
              identical(names(fwrcdr_widths),
                        names(IMGT_DEFAULT_FWRCDR_WIDTHS)),
              is.character(region_types))

    extra_cols <- mcols(gapped_V_alleles, use.names=FALSE)
    if (length(extra_cols) != 0L)
        mcols(gapped_V_alleles) <- NULL

    fwrcdr_ends <- cumsum(fwrcdr_widths)
    end <- pmin(fwrcdr_ends[["fwr1"]], width(gapped_V_alleles))
    fwr1 <- subseq(gapped_V_alleles, end=end)
    start <- end + 1L
    end <- pmin(fwrcdr_ends[["cdr1"]], width(gapped_V_alleles))
    cdr1 <- subseq(gapped_V_alleles, start=start, end=end)
    start <- end + 1L
    end <- pmin(fwrcdr_ends[["fwr2"]], width(gapped_V_alleles))
    fwr2 <- subseq(gapped_V_alleles, start=start, end=end)
    start <- end + 1L
    end <- pmin(fwrcdr_ends[["cdr2"]], width(gapped_V_alleles))
    cdr2 <- subseq(gapped_V_alleles, start=start, end=end)
    start <- end + 1L
    end <- pmin(fwrcdr_ends[["fwr3"]], width(gapped_V_alleles))
    fwr3 <- subseq(gapped_V_alleles, start=start, end=end)
    start <- end + 1L
    cdr3 <- subseq(gapped_V_alleles, start=start)

    ans <- DataFrame(fwr1=fwr1, cdr1=cdr1,
                     fwr2=fwr2, cdr2=cdr2,
                     fwr3=fwr3, cdr3=cdr3)
    ans <- ans[ , region_types, drop=FALSE]
    ans <- cbind(allele_name=names(gapped_V_alleles), ans)
    if (length(extra_cols) != 0L)
        ans <- cbind(ans, extra_cols)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .print_Vparts()
###

.extract_Vparts_extra_cols <- function(Vparts, region_types)
{
    stopifnot(is(Vparts, "DataFrame"), is.character(region_types))
    core_colnames <- c("allele_name", region_types)
    core_col_idx <- match(core_colnames, colnames(Vparts))
    stopifnot(!anyNA(core_col_idx))
    Vparts[-core_col_idx]
}

.make_Vparts_header_line <- function(labelW, part_widths, sep, extra_colnames)
{
    stopifnot(isSingleInteger(labelW),
              is.integer(part_widths),
              all(names(part_widths) %in% .VPART_NAMES),
              isSingleString(sep),
              is.null(extra_colnames) || is.character(extra_colnames))
    arrows <- vapply(seq_along(part_widths),
        function(i) {
            text <- names(part_widths)[[i]]
            if (i == length(part_widths)) {
                Rarrow_char <- "-"
                justify <- "center-right"
            } else {
                Rarrow_char <- ">"
                justify <- "center-left"
            }
            make_ascii_double_arrow(part_widths[[i]], text,
                                    justify=justify, Rarrow_char=Rarrow_char)
        },
        character(1))
    ans <- paste0(strrep(" ", labelW), paste(arrows, collapse=sep))
    if (length(extra_colnames) != 0L)
        ans <- paste0(ans, " ", paste(extra_colnames, collapse=" "))
    ans
}

.make_Vparts_line <- function(i, labels, Vparts, Rfillers, sep, extra_cols)
{
    stopifnot(isSingleInteger(i),
              is.character(labels),
              is(Vparts, "DataFrame"),
              is.matrix(Rfillers),
              all(colnames(Rfillers) %in% .VPART_NAMES),
              isSingleString(sep),
              is.matrix(extra_cols),
              length(labels) == nrow(Vparts),
              length(labels) == nrow(Rfillers),
              length(labels) == nrow(extra_cols))
    colored_strings <- vapply(colnames(Rfillers),
        function(region_type)
            paste0(paint_XString(Vparts[ , region_type][[i]]),
                   Rfillers[i, region_type]),
        character(1))
    ans <- paste0(labels[[i]], paste(colored_strings, collapse=sep))
    if (ncol(extra_cols) != 0L)
        ans <- paste0(ans, " ", paste(extra_cols[i, ], collapse=" "))
    ans
}

.print_Vparts <- function(Vparts, fwrcdr_widths, region_types,
                          filler=".", sep=" ")
{
    extra_cols <- .extract_Vparts_extra_cols(Vparts, region_types)
    stopifnot(is.integer(fwrcdr_widths),
              identical(names(fwrcdr_widths),
                        names(IMGT_DEFAULT_FWRCDR_WIDTHS)))

    ## Prepare 'labels'.
    allele_names <- Vparts[ , "allele_name"]
    labels <- paste0(format(seq_len(nrow(Vparts))), ". ",
                     format(paste0(allele_names, ": ")))
    labelW <- nchar(labels[[1L]])
    stopifnot(all(nchar(labels) == labelW))

    ## Complete 'part_widths'.
    if ("cdr3" %in% region_types)
        fwrcdr_widths <- c(fwrcdr_widths, cdr3=max(width(Vparts[ , "cdr3"])))
    part_widths <- fwrcdr_widths[region_types]

    ## Compute 'Rfillers'.
    Rfillers <- lapply(seq_along(part_widths),
        function(i) {
            region_type <- region_types[[i]]
            strrep(filler, part_widths[[i]] - width(Vparts[ , region_type]))
        })
    Rfillers <- do.call(cbind, setNames(Rfillers, region_types))

    ## Format 'extra_cols'.
    extra_cols <- DataFrame_as_formatted_matrix(extra_cols)
    extra_colnames <- colnames(extra_cols)

    ## Print the header.
    line <- .make_Vparts_header_line(labelW, part_widths, sep, extra_colnames)
    stopifnot(isSingleString(line))
    message(line)

    ## Print the alleles sequences.
    locus <- substr(allele_names, 1L, 3L)
    group_lens <- runLength(Rle(locus))
    if (length(group_lens) == length(unique(locus))) {
        groupend_idx <- cumsum(head(group_lens, n=-1L))
    } else {
        groupend_idx <- integer(0)
    }
    for (i in seq_len(nrow(Vparts))) {
        line <- .make_Vparts_line(i, labels, Vparts, Rfillers, sep, extra_cols)
        stopifnot(isSingleString(line))
        message(line)
        if (i %in% groupend_idx)
            message("")
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### print_gapped_V_alleles()
###

print_gapped_V_alleles <- function(gapped_V_alleles,
                                   fwrcdr_widths=IMGT_DEFAULT_FWRCDR_WIDTHS,
                                   translate=FALSE,
                                   region_types=NULL,
                                   igblast_organism=NA, filler=".", sep=" ")
{
    ## normarg_gapped_V_alleles() can return a DNAStringSet or BStringSet
    ## object.
    gapped_V_alleles <- normarg_gapped_V_alleles(gapped_V_alleles)
    if (!is(gapped_V_alleles, "DNAStringSet"))
        gapped_V_alleles <- as(gapped_V_alleles, "DNAStringSet")
    if (is.matrix(fwrcdr_widths))
        stop(wmsg("'fwrcdr_widths' cannot be a matrix"))
    fwrcdr_widths <- normarg_fwrcdr_widths(fwrcdr_widths)
    if (!isTRUEorFALSE(translate))
        stop(wmsg("'translate' must be TRUE or FALSE"))
    region_types <- .normarg_region_types2(region_types)
    if (!isSingleString(filler) || nchar(filler) != 1L)
        stop(wmsg("'filler' must be a single character"))
    if (!isSingleString(sep))
        stop(wmsg("'sep' must be a single string"))

    if (translate) {
        gapped_V_alleles <- translate_codons(gapped_V_alleles)
    } else {
        fwrcdr_widths <- fwrcdr_widths * 3L
    }
    Vparts <- .extract_Vparts(gapped_V_alleles, fwrcdr_widths, region_types,
                              igblast_organism=igblast_organism)
    .print_Vparts(Vparts, fwrcdr_widths, region_types, filler=filler, sep=sep)
    invisible(Vparts)
}


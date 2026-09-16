
.get_all_fwrcdr_widths <- function(loci_prefix)
{
    imgt_additional_insertions <- switch(loci_prefix,
        IG=igblastr:::.IMGT_IG_ADDITIONAL_INSERTIONS,
        TR=igblastr:::.IMGT_TR_ADDITIONAL_INSERTIONS,
        stop(wmsg("invalid 'loci_prefix'")))
    imgt_organisms <- names(imgt_additional_insertions)
    all_fwrcdr_widths <- lapply(setNames(imgt_organisms, imgt_organisms),
                                get_fwrcdr_widths_for_imgt_organism,
                                loci_prefix=loci_prefix)
    all_storage_modes <- vapply(all_fwrcdr_widths, storage.mode, character(1))
    expect_true(all(all_storage_modes == "integer"))
    all_fwrcdr_widths
}

test_that("get_fwrcdr_widths_for_imgt_organism()", {

    fwrdcr_names <- names(IMGT_DEFAULT_FWRCDR_WIDTHS)

    ## --- 'loci_prefix' set to "IG" ---

    all_fwrcdr_widths <- .get_all_fwrcdr_widths("IG")
    all_classes <- vapply(all_fwrcdr_widths,
        function(fwrcdr_widths) class(fwrcdr_widths)[[1L]],
        character(1)
    )
    expect_true(all(all_classes %in% c("integer", "matrix")))

    ## Check list elements of class "integer".
    all_fwrcdr_widths1 <- all_fwrcdr_widths[all_classes == "integer"]
    expect_identical(names(all_fwrcdr_widths1),
                     c("Bos_taurus",
                       "Gorilla_gorilla_gorilla",
                       "Homo_sapiens",
                       "Lemur_catta",
                       "Mus_musculus",
                       "Oryctolagus_cuniculus",
                       "Sus_scrofa"))
    ok <- vapply(all_fwrcdr_widths1,
        function(fwrcdr_widths)
            identical(names(fwrcdr_widths), fwrdcr_names),
        logical(1))
    expect_true(all(ok))

    ## Check list elements of class "matrix".
    all_fwrcdr_widths2 <- all_fwrcdr_widths[all_classes == "matrix"]
    ok <- vapply(all_fwrcdr_widths2,
        function(fwrcdr_widths)
            identical(colnames(fwrcdr_widths), fwrdcr_names),
        logical(1))
    expect_true(all(ok))
    valid_groups <- paste0(igblastr:::IG_LOCI, "V")
    ok <- vapply(all_fwrcdr_widths2,
        function(fwrcdr_widths) {
            m <- match(rownames(fwrcdr_widths), valid_groups)
            isFALSE(is.unsorted(m, strict=TRUE))
            }, logical(1))
    expect_true(all(ok))

    ## More checks.
    expect_true(is.na(get_fwrcdr_widths_for_imgt_organism("Sasquatch")))
    fwrcdr_widths <- get_fwrcdr_widths_for_imgt_organism("Homo_sapiens",
                                                         loci_prefix="IG")
    expect_identical(fwrcdr_widths, IMGT_DEFAULT_FWRCDR_WIDTHS)
    fwrcdr_widths <- get_fwrcdr_widths_for_imgt_organism("Gallus_gallus",
                                                         loci_prefix="IG")
    expected <- rbind(IGHV=IMGT_DEFAULT_FWRCDR_WIDTHS,
                      IGLV=c(27L, 12L, 18L, 10L, 39L))
    expect_identical(fwrcdr_widths, expected)

    ## --- 'loci_prefix' set to "TR" ---

    all_fwrcdr_widths <- .get_all_fwrcdr_widths("TR")
    all_classes <- vapply(all_fwrcdr_widths,
        function(fwrcdr_widths) class(fwrcdr_widths)[[1L]],
        character(1)
    )
    expect_true(all(all_classes %in% c("integer", "matrix")))

    ## Check list elements of class "integer".
    all_fwrcdr_widths1 <- all_fwrcdr_widths[all_classes == "integer"]
    expect_identical(names(all_fwrcdr_widths1),
                     c("Camelus_dromedarius",
                       "Homo_sapiens",
                       "Oryctolagus_cuniculus",
                       "Pongo_abelii"))
    ok <- vapply(all_fwrcdr_widths1,
        function(fwrcdr_widths)
            identical(names(fwrcdr_widths), fwrdcr_names),
        logical(1))
    expect_true(all(ok))

    ## Check list elements of class "matrix".
    all_fwrcdr_widths2 <- all_fwrcdr_widths[all_classes == "matrix"]
    ok <- vapply(all_fwrcdr_widths2,
        function(fwrcdr_widths)
            identical(colnames(fwrcdr_widths), fwrdcr_names),
        logical(1))
    expect_true(all(ok))
    valid_groups <- paste0(igblastr:::TR_LOCI, "V")
    ok <- vapply(all_fwrcdr_widths2,
        function(fwrcdr_widths) {
            m <- match(rownames(fwrcdr_widths), valid_groups)
            isFALSE(is.unsorted(m, strict=TRUE))
            }, logical(1))
    expect_true(all(ok))

    ## More checks.
    expect_true(is.na(get_fwrcdr_widths_for_imgt_organism("Sasquatch",
                                                          loci_prefix="TR")))
    fwrcdr_widths <- get_fwrcdr_widths_for_imgt_organism("Pongo_abelii",
                                                          loci_prefix="TR")
    expect_identical(fwrcdr_widths, IMGT_DEFAULT_FWRCDR_WIDTHS)
    fwrcdr_widths <- get_fwrcdr_widths_for_imgt_organism("Danio_rerio",
                                                         loci_prefix="TR")
    expected <- rbind(TRAV=c(fwr1=27L, cdr1=12L, fwr2=18L, cdr2=12L, fwr3=41L))
    expect_identical(fwrcdr_widths, expected)
})


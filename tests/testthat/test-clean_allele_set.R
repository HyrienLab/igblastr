test_that(".drop_repeated_alleles()", {
    drop_repeated_alleles <- igblastr:::.drop_repeated_alleles

    x <- c(X="ACTT", Y="ACTT", X="CC", Y="TT", Y="GGGGGGGG",
           X="ACTT", X="ACTT", Z="", Y="TT", Y="GgGgGgGg")

    expect_identical(drop_repeated_alleles(x), x[-c(6:7, 9)])

    allele_set <- DNAStringSet(x)
    object <- drop_repeated_alleles(allele_set)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), x[-c(6:7, 9:10)])

    y <- rev(x)
    expect_identical(drop_repeated_alleles(y), y[-c(5L, 7L, 10L)])

    allele_set <- DNAStringSet(y)
    object <- drop_repeated_alleles(allele_set)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), toupper(y[-c(5:7, 10L)]))
})

.check_object_mcols <- function(object)
{
    object_mcols <- mcols(object)
    expect_true(is(object_mcols, "DataFrame"))
    expected_colnames <- c("locus",
                           igblastr:::EXTENDED_V_GENE_DELINEATION_COLNAMES,
                           "chain_type")
    expect_identical(colnames(object_mcols), expected_colnames)
    expect_identical(object_mcols[ , "allele_name"], names(object))
}

test_that("clean_V_allele_set()", {
    clean_V_allele_set <- igblastr:::clean_V_allele_set

    ## ======== WITHOUT GAPS ========

    allele_set0 <- DNAStringSet(c(A1="GG", A2="AAAAACG", A3="T",
                                  A4="GG", A5="CTAATA",  A6="GG",
                                  A7="TTT", A8="GG", A9="TTTTTTTTTT"))
    object <- clean_V_allele_set(allele_set0)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(allele_set0))

    allele_set <- allele_set0
    names(allele_set) <- paste0("   ", names(allele_set), "   \t ")
    object <- clean_V_allele_set(allele_set)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(allele_set0))

    allele_set <- allele_set0
    names(allele_set) <- paste0("\t stuff |",
                                names(allele_set),
                                "\t   |more|| stuff  ")
    object <- clean_V_allele_set(allele_set)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(allele_set0))

    allele_set <- allele_set0
    names(allele_set)[2L] <- NA
    expect_error(clean_V_allele_set(allele_set), regexp="NAs")

    allele_set <- allele_set0
    names(allele_set)[2L] <- "   \t\t"
    expect_error(clean_V_allele_set(allele_set), regexp="empty")

    allele_set <- allele_set0
    names(allele_set)[2L] <- "Z"
    expect_error(clean_V_allele_set(allele_set),
                 regexp="less than 2-character long")

    ## With "repeated" alleles (i.e. alleles with identical DNA
    ## sequences **and** names) but no "ambiguous" alleles (i.e.
    ## alleles with same name but different sequences).
    allele_set <- allele_set0
    names(allele_set)[8L] <- "A6"
    object <- clean_V_allele_set(allele_set)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(allele_set[-8L]))

    ## With "repeated" alleles (i.e. alleles with identical DNA
    ## sequences **and** names) AND "ambiguous" alleles (i.e.
    ## alleles with same name but different sequences).
    allele_set <- allele_set0
    names(allele_set)[c(1:2, 9L)] <- "A1"
    names(allele_set)[c(3L, 5:6, 8L)] <- "A3"
    expected_names <- c("A1a", "A1b", "A3a", "A4", "A3b", "A3c", "A7", "A1c")
    expected_object <- setNames(allele_set0[-8L], expected_names)
    errmsg <- "The following allele names are ambiguous: A1, A3"
    expect_error2(clean_V_allele_set(allele_set), errmsg)
    object <- clean_V_allele_set(allele_set, disambiguate.allele.names=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(expected_object))
    expect_true(is.null(mcols(object)))

    regexp <- "V allele sequences have no gaps"
    expect_warning(
        object <- clean_V_allele_set(allele_set, gapped=TRUE,
                                     disambiguate.allele.names=TRUE),
        regexp=regexp
    )
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(expected_object))
    expect_true(is.null(mcols(object)))

    ## ======== WITH GAPS ========

    allele_set1 <- DNAStringSet(c(A1="AC.G..T", A2="AAA...AACGT", A3="G...G"))
    errmsg <- "Some V allele sequences have gaps"
    expect_error2(clean_V_allele_set(allele_set1), errmsg)

    object <- clean_V_allele_set(allele_set1, gapped=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), c(A1="ACGT", A2="AAAAACGT", A3="GG"))
    expect_true(is.null(mcols(object)))

    errmsg <- paste0("'with.intdata=TRUE' can only be used ",
                     "when 'gapped' is TRUE")
    expect_error2(clean_V_allele_set(allele_set1, with.intdata=TRUE), errmsg)

    errmsg <- paste0("'allele_set' must have a \"locus\" ",
                     "metadata column when 'with.intdata' is TRUE")
    expect_error2(
        clean_V_allele_set(allele_set1, gapped=TRUE, with.intdata=TRUE),
        errmsg
    )

    mcols(allele_set1)$locus <- "IGH"
    allele_set <- allele_set1
    names(allele_set) <- paste0("\t stuff |",
                                names(allele_set),
                                "\t   |more|| stuff  ")
    object <- clean_V_allele_set(allele_set, gapped=TRUE, with.intdata=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), c(A1="ACGT", A2="AAAAACGT", A3="GG"))
    .check_object_mcols(object)

    allele_set2 <- DNAStringSet(c(A1="ACGT", A2="AAAAACGT", A3="G...G"))
    errmsg <- "Some V allele sequences have gaps"
    expect_error2(clean_V_allele_set(allele_set2), errmsg)

    allele_set <- allele_set2
    names(allele_set) <- paste0("\t ", names(allele_set), "\t    ")
    regexp <- "V allele sequences have no gaps"
    expect_warning(
        object <- clean_V_allele_set(allele_set, gapped=TRUE),
        regexp=regexp
    )
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), c(A1="ACGT", A2="AAAAACGT", A3="GG"))
    expect_true(is.null(mcols(object)))

    mcols(allele_set2)$locus <- "IGK"
    allele_set <- allele_set2
    names(allele_set) <- paste0("\t ", names(allele_set), "\t    ")
    regexp <- "V allele sequences have no gaps"
    expect_warning(
        object <- clean_V_allele_set(allele_set,
                                     gapped=TRUE, with.intdata=TRUE),
        regexp=regexp
    )
    expect_true(is(object, "DNAStringSet"))
    .check_object_mcols(object)

    allele_set3 <- DNAStringSet(c(A1="G.G.", A2="AAA........AACG",
                                  A3="T..", A4="....GG..", A5="CT..AATA",
                                  A6="G....G", A7="T..T..T", A8=".G...G",
                                  A9=".........TTT...T.T.T.T..TT.....T"))
    mcols(allele_set3)$locus <- "TRA"

    object <- clean_V_allele_set(allele_set3, gapped=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(allele_set0))
    expect_identical(colnames(mcols(object)), colnames(mcols(allele_set3)))

    object <- clean_V_allele_set(allele_set3, gapped=TRUE, with.intdata=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(allele_set0))
    .check_object_mcols(object)

    ## With "repeated" alleles (i.e. alleles with identical **ungapped** DNA
    ## sequences **and** names) AND "ambiguous" alleles (i.e.
    ## alleles with same name but different **ungapped** sequences).
    allele_set <- allele_set3
    names(allele_set)[c(1:2, 9L)] <- "A1"
    names(allele_set)[c(3L, 5:6, 8L)] <- "A3"
    expected_names <- c("A1a", "A1b", "A3a", "A4", "A3b", "A3c", "A7", "A1c")
    expected_object <- setNames(allele_set0[-8L], expected_names)
    object <- clean_V_allele_set(allele_set, gapped=TRUE,
                                 disambiguate.allele.names=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(expected_object))
    expect_identical(colnames(mcols(object)), colnames(mcols(allele_set3)))

    ## Fails because the gaps in A6 and A8 imply different coding_frame_start.
    errmsg <- paste0("alleles with identical ungapped sequences ",
                     "and names must have gaps that result in ",
                     "identical annotations")
    expect_error2(
        clean_V_allele_set(allele_set, gapped=TRUE, with.intdata=TRUE),
        errmsg
    )

    ## Let's put gaps in A6 and A8 that result in the two alleles having
    ## identical annotations.
    allele_set <- allele_set3
    allele_set[c("A6", "A8")] <- c(".......GG", ".G....G..")
    mcols(allele_set)$locus <- "TRB"
    names(allele_set)[c(1:2, 9L)] <- "A1"
    names(allele_set)[c(3L, 5:6, 8L)] <- "A3"
    object <- clean_V_allele_set(allele_set, gapped=TRUE, with.intdata=TRUE,
                                 disambiguate.allele.names=TRUE)
    expect_identical(as.character(object), as.character(expected_object))
    .check_object_mcols(object)
})


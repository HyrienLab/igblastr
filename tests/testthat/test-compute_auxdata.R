
test_that("compute_auxdata()", {

    ## Un unrealistic J allele sequence that contains 1 occurrence of
    ## WGXG and 2 occurrences of FGXG.
    dna <- "TACAGCAACTGGGGACTTGGTTTGGGGAAGGGATTCGGAGTTGGCGTGTTACCAA"
    J_alleles <- DNAStringSet(c(`IGHJ5*01`=dna, `IGKJ2*03`=dna))

    computed_auxdata <- compute_auxdata(J_alleles)
    igblastr:::check_auxdata_col2class(computed_auxdata)
    expected <- data.frame(
        allele_name=names(J_alleles),
        coding_frame_start=c(0L, 2L),
        chain_type=c("JH", "JK"),
        cdr3_end=c(8L, 19L),
        extra_bps=c(1L, 2L)
    )
    expect_identical(computed_auxdata, expected)
    codon_starts <- c(`IGHJ5*01`=1, `IGKJ2*03`=3)
    expect_identical(compute_auxdata(J_alleles, codon_starts=codon_starts),
                     expected)

    codon_starts <- c(`IGHJ5*01`=2, `IGKJ2*03`=1)
    expected <- data.frame(
        allele_name=names(J_alleles),
        coding_frame_start=c(1L, 0L),
        chain_type=c("JH", "JK"),
        cdr3_end=c(NA_integer_, 32L),
        extra_bps=c(0L, 1L)
    )
    expect_identical(compute_auxdata(J_alleles, codon_starts=codon_starts),
                     expected)
    codon_starts <- c(`IGHJ5*01`=3, `IGKJ2*03`=2)
    expected <- data.frame(
        allele_name=names(J_alleles),
        coding_frame_start=c(2L, 1L),
        chain_type=c("JH", "JK"),
        cdr3_end=c(NA_integer_, NA_integer_),
        extra_bps=c(2L, 0L)
    )
    expect_identical(compute_auxdata(J_alleles, codon_starts=codon_starts),
                     expected)

    ## --- for human J alleles from OGRDB ---

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
    J_alleles <- load_germline_sequences(db_name, region_types="J")
    computed_auxdata <- compute_auxdata(J_alleles)
    igblastr:::check_auxdata_col2class(computed_auxdata)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Now we're going to check that 'computed_auxdata' agrees with the
    ## auxiliary data included in IgBLAST. More precisely, we're going to
    ## check that it's a subset of 'load_auxdata("human", which="original")'.

    ## All the J alleles in _OGRDB.human.IGH+IGK+IGL.202605 are annotated
    ## in human_gl.aux so we expect no NAs in 'm' below.
    orig_auxdata <- igblastr:::load_and_fix_igblast_auxdata("human")
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    expect_false(anyNA(m))
    orig_auxdata <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m)
    expect_identical(computed_auxdata, orig_auxdata)
})

test_that("find_discordant_auxdata()/complete_auxdata_with_ref()", {
    find_discordant_auxdata <- igblastr:::find_discordant_auxdata
    complete_auxdata_with_ref <- igblastr:::complete_auxdata_with_ref

    ## --- trivial case ---

    auxdata1 <- data.frame(
        allele_name=c("A", "B", "C", "D", "E", "F"),
        coding_frame_start=c(1L, 2L, 1L, 0L, 1L, 1L),
        chain_type=c("JH", "JH", "JK", "JL", "JH", "JK"),
        cdr3_end=c(11:13, NA, 15:16),
        extra_bps=1L
    )
    expect_identical(complete_auxdata_with_ref(auxdata1, auxdata1),
                     auxdata1)

    ## --- discordant data ---

    auxdata2 <- data.frame(
        allele_name=c("E", "F", "G", "H", "B", "D", "I", "A"),
        coding_frame_start=c(NA, 1L, 1L, 2L, 1L, 0L, 0L, 1L),
        chain_type=c("JH", "JK", "JH", "JH", "JH", "JL", "JL", "JH"),
        cdr3_end=c(11L, 11:13, 12L, 15:16, 11L),
        extra_bps=1L
    )

    expect_error(complete_auxdata_with_ref(auxdata1, auxdata2),
                 regexp="discordant")
    disc_rowpairs <- find_discordant_auxdata(auxdata1, auxdata2)
    expected <- data.frame(rowidx1=c(2L, 5L, 6L), rowidx2=c(5L, 1L, 2L))
    expect_identical(disc_rowpairs, expected)

    nopairs <- data.frame(rowidx1=integer(0), rowidx2=integer(0))
    disc_rowpairs <- find_discordant_auxdata(auxdata1,
                                             auxdata2[-disc_rowpairs[[2L]], ])
    expect_identical(disc_rowpairs, nopairs)
    disc_rowpairs <- find_discordant_auxdata(auxdata1[-disc_rowpairs[[1L]], ],
                                             auxdata2)
    expect_identical(disc_rowpairs, nopairs)

    ## Swapping the two data.frames returns the same pairs but they are
    ## possibly in a different order.
    expect_error(complete_auxdata_with_ref(auxdata2, auxdata1),
                 regexp="discordant")
    disc_rowpairs <- find_discordant_auxdata(auxdata2, auxdata1)
    expected <- data.frame(rowidx1=c(1L, 2L, 5L), rowidx2=c(5L, 6L, 2L))
    expect_identical(disc_rowpairs, expected)

    ## --- concordant data ---

    auxdata1[c(2L, 5L), "coding_frame_start"] <- NA
    auxdata2[1:2, "cdr3_end"] <- NA

    disc_rowpairs <- find_discordant_auxdata(auxdata1, auxdata2)
    expect_identical(disc_rowpairs, nopairs)
    expected <- auxdata1
    expected[2L, "coding_frame_start"] <- 1L
    expected[4L, "cdr3_end"] <- 15L
    auxdata1a <- complete_auxdata_with_ref(auxdata1, auxdata2)
    expect_identical(auxdata1a, expected)
    expect_identical(complete_auxdata_with_ref(auxdata1a, auxdata2),
                     auxdata1a)

    disc_rowpairs <- find_discordant_auxdata(auxdata2, auxdata1)
    expect_identical(disc_rowpairs, nopairs)
    expected <- auxdata2
    expected[c(1L, 2L), "cdr3_end"] <- c(15L, 16L)
    auxdata2a <- complete_auxdata_with_ref(auxdata2, auxdata1)
    expect_identical(auxdata2a, expected)
    expect_identical(complete_auxdata_with_ref(auxdata2a, auxdata1),
                     auxdata2a)
})

test_that(".get_fwr4_start_from_fwr4refset_aa_dist", {
    get_fwr4_start_from_fwr4refset_aa_dist <-
        igblastr:::.get_fwr4_start_from_fwr4refset_aa_dist

    fwr4refset <- AAStringSet(c(
        "WGTGTLVTISS",
        "FGAGTKVEIK",
        "FGPGTKLIITS",
        "IVVSPKVEVVIG",
        "FGSRTKVIVT",
        "FKGGFGTLATIS",
        "FAEGTVLVVTSP",
        "FEKGTYLEEQQ",
        "FGIGTKLQVIP",
        "LVTIGGTKVS",
        "TGTLVTISISS"
    ))

    ## Min distance is 0 (exact match) and is achieved for P = 8.
    ## Second best distance is 6 and is achieved for P = 6.
    aa_string <- AAString("LGKKIKVFGPGTKLIIT")
    stopifnot(subseq(aa_string, start=8L, width=10L) ==
              subseq(fwr4refset[[3L]], start=1L, width=10L))
    P <- get_fwr4_start_from_fwr4refset_aa_dist(aa_string, fwr4refset)
    expect_identical(P, 8L)
    P <- get_fwr4_start_from_fwr4refset_aa_dist(aa_string, fwr4refset,
                                                standout.by=6)
    expect_identical(P, 8L)
    P <- get_fwr4_start_from_fwr4refset_aa_dist(aa_string, fwr4refset,
                                                standout.by=7)
    expect_identical(P, NA_integer_)

    ## Min distance is 3 and is achieved for P = 5.
    ## Second best distance is 6 and is achieved for P = 7 and P = 12.
    aa_string <- AAString("QQQQSAEGTKLIVTSVVVPPTVGG")
    P <- get_fwr4_start_from_fwr4refset_aa_dist(aa_string, fwr4refset)
    expect_identical(P, NA_integer_)
    P <- get_fwr4_start_from_fwr4refset_aa_dist(aa_string, fwr4refset,
                                                max.dist=3)
    expect_identical(P, NA_integer_)
    P <- get_fwr4_start_from_fwr4refset_aa_dist(aa_string, fwr4refset,
                                                max.dist=3, standout.by=3)
    expect_identical(P, 5L)
    P <- get_fwr4_start_from_fwr4refset_aa_dist(aa_string, fwr4refset,
                                                max.dist=5, standout.by=3)
    expect_identical(P, 5L)
})

test_that(".get_fwr4_start_from_fwr4refset_dna_dist", {
    get_fwr4_start_from_fwr4refset_dna_dist <-
        igblastr:::.get_fwr4_start_from_fwr4refset_dna_dist

    fwr4refset <- DNAStringSet(c(
        "TTACACCTGGAGATTAGGGAAAGTAGGGTC",     # 30 nucleotides
        "ACTAATTAAAGGAGCTTATTTAAAAGA",        # 27 nucleotides
        "GCTAGGATGGATTGATGTGTAGTAGGTCAAAT",   # 32 nucleotides
        "AGCTAGGATGGACTGATCTGTAGTAGGTCA",     # 30 nucleotides
        "AGGATGGAcTGATcTGTAGTAGGTCAATTTgggg"  # 34 nucleotides
    ))

    ## Min distance is 2 and is achieved for P = 10.
    ## The exact match at P = 9  is ignored because we only look at
    ## positions 1, 4, 7, 10, etc..
    dna_string <- DNAString("aaagggaaaGCTAGGATGGAcTGATcTGTAGTAGGTCAA")
    stopifnot(subseq(dna_string, start=9L, width=30L) == fwr4refset[[4L]])
    P <- get_fwr4_start_from_fwr4refset_dna_dist(dna_string, fwr4refset)
    expect_identical(P, 10L)
    P <- get_fwr4_start_from_fwr4refset_dna_dist(dna_string, fwr4refset,
                                                 max.dist=2)
    expect_identical(P, 10L)
    P <- get_fwr4_start_from_fwr4refset_dna_dist(dna_string, fwr4refset,
                                                 max.dist=1)
    expect_identical(P, NA_integer_)

    ## Min distance is 0 (exact match) and is achieved for P = 13.
    ## Second best distance is 5 and is achieved for P = 10.
    dna_string <- DNAString("aaagggcccaaaAGGATGGAcTGATcTGTAGTAGGTCAATTTcgcg")
    stopifnot(subseq(dna_string, start=13L, width=30L) ==
              subseq(fwr4refset[[5L]], start=1L, width=30L))
    P <- get_fwr4_start_from_fwr4refset_dna_dist(dna_string, fwr4refset,
                                                 max.dist=0, standout.by=5)
    expect_identical(P, 13L)
    P <- get_fwr4_start_from_fwr4refset_dna_dist(dna_string, fwr4refset,
                                                 max.dist=0, standout.by=6)
    expect_identical(P, NA_integer_)

    ## Min distance is 3 and is achieved for P = 7.
    dna_string <- DNAString("aaaaaaACTAATTAAAGGAGCTTATTTAAAAGAaaa")
    stopifnot(subseq(dna_string, start=7L, width=27L) == fwr4refset[[2L]])
    P <- get_fwr4_start_from_fwr4refset_dna_dist(dna_string, fwr4refset)
    expect_identical(P, 7L)
    P <- get_fwr4_start_from_fwr4refset_dna_dist(dna_string, fwr4refset,
                                                 max.dist=3)
    expect_identical(P, 7L)
    P <- get_fwr4_start_from_fwr4refset_dna_dist(dna_string, fwr4refset,
                                                 max.dist=2)
    expect_identical(P, NA_integer_)
})

test_that(".get_fwr4_start_from_fwr4refset_dna_PWM", {
    get_fwr4_start_from_fwr4refset_dna_PWM <-
        igblastr:::.get_fwr4_start_from_fwr4refset_dna_PWM

    fwr4refset <- DNAStringSet(c(
        "TTACACCTGGAGATTAGGGAAAGTAGGGTC",     # 30 nucleotides
        "GCTAGGATGGATTGATGTGTAGTAGGTCAAAT",   # 32 nucleotides
        "AGCTAGGATGGACTGATCTGTAGTAGGTCA",     # 30 nucleotides
        "AGGATGGAcTGATcTGTAGTAGGTCAATTTgggg"  # 34 nucleotides
    ))
    pwm <- PWM(subseq(fwr4refset, start=1L, end=30L))

    dna_string <- DNAString("aaagggaaaGCTAGGATGGAcTGATcTGTAGTAGGTCAA")
    P <- get_fwr4_start_from_fwr4refset_dna_PWM(dna_string, pwm)
    expect_identical(P, NA_integer_)
    P <- get_fwr4_start_from_fwr4refset_dna_PWM(dna_string, pwm,
                                                standout.by=0.17)
    expect_identical(P, 10L)
})

test_that("solve_cdr3_ends_using_fwr4_*_comparisons", {
    ## The auxdata for IMGT rat has 1 unsolved allele.
    db_name <- install_IMGT_germline_db("202631-1", "Rattus norvegicus",
                                        overwrite=TRUE)
    auxdata0 <- load_auxdata(db_name)
    solve_me <- is.na(auxdata0[ , "cdr3_end"])
    stopifnot(identical(auxdata0[solve_me, "allele_name"], "IGKJ3*01"))
    J_alleles <- load_germline_sequences(db_name, region_types="J")

    auxdata <- solve_cdr3_ends_using_fwr4_aa_comparisons(auxdata0, J_alleles)
    expect_identical(auxdata, auxdata0)
    auxdata <- solve_cdr3_ends_using_fwr4_dna_comparisons(auxdata0, J_alleles)
    expect_identical(auxdata, auxdata0)
    auxdata <- solve_cdr3_ends_using_fwr4_dna_PWM(auxdata0, J_alleles)
    expect_identical(auxdata, auxdata0)

    auxdata <- solve_cdr3_ends_using_fwr4_aa_comparisons(auxdata0, J_alleles,
                                                         max.dist=5)
    expect_identical(auxdata, auxdata0)
    auxdata <- solve_cdr3_ends_using_fwr4_dna_comparisons(auxdata0, J_alleles,
                                                          max.dist=9)
    expect_identical(auxdata, auxdata0)
    auxdata <- solve_cdr3_ends_using_fwr4_dna_PWM(auxdata0, J_alleles,
                                                  min.score=0.82,
                                                  standout.by=0.37)
    expect_identical(auxdata, auxdata0)
    auxdata <- solve_cdr3_ends_using_fwr4_dna_PWM(auxdata0, J_alleles,
                                                  min.score=0.81,
                                                  standout.by=0.38)
    expect_identical(auxdata, auxdata0)

    auxdata1 <- solve_cdr3_ends_using_fwr4_aa_comparisons(auxdata0, J_alleles,
                                                          max.dist=5,
                                                          standout.by=4)
    expect_identical(auxdata1[solve_me , "cdr3_end"], 6L)
    auxdata2 <- solve_cdr3_ends_using_fwr4_dna_comparisons(auxdata0, J_alleles,
                                                           max.dist=9,
                                                           standout.by=11)
    expect_identical(auxdata2, auxdata1)
    auxdata3 <- solve_cdr3_ends_using_fwr4_dna_PWM(auxdata0, J_alleles,
                                                   min.score=0.81,
                                                   standout.by=0.37)
    expect_identical(auxdata3, auxdata1)

    DF <- suppressMessages(print_J_alleles(J_alleles, auxdata1, translate=TRUE))
    expect_true(DF[solve_me , "fwr4"] == AAStringSet("ISDETRLEIK"))
})


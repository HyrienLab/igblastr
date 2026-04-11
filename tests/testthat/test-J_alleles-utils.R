### Install IMGT germline dbs used in tests below. Note that only the first
### installation actually triggers a download from IMGT. All subsequent
### installations obtain the data from the IMGT local store (located
### in 'igblastr_cache(IMGT_STORE)') so are very fast and work offline.
install_IMGT_germline_db("202614-2", "Homo sapiens",
                         without.auxdata=TRUE, overwrite=TRUE)
install_IMGT_germline_db("202614-2", "Mus musculus",
                         without.auxdata=TRUE, overwrite=TRUE)
install_IMGT_germline_db("202614-2", "Rattus norvegicus",
                         without.auxdata=TRUE, overwrite=TRUE)
install_IMGT_germline_db("202614-2", "Oryctolagus cuniculus",
                         without.auxdata=TRUE, overwrite=TRUE)

test_that("translate_J_alleles()", {

    ## --- for human J alleles (from OGRDB and IMGT) ---

    auxdata <- load_and_fix_human_auxdata()

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    J_aa <- translate_J_alleles(J_alleles, auxdata)
    expect_true(is.character(J_aa))
    expect_identical(names(J_aa), names(J_alleles))
    expect_false(anyNA(J_aa))
    expect_identical(J_aa[["IGHJ1*01"]], "AEYFQHWGQGTLVTVSS")
    expect_true(all(grepl("[WF]G.G", J_aa)))

    db_name <- "IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    ## 4 human J alleles in IMGT release 202614-2 have no entries
    ## in 'auxdata'. Note that this could change in the future.
    expected_unknown <- c("IGHJ5*03", "IGHJ5*04", "IGKJ4*03", "IGLJ2A*01")
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    expect_identical(names(J_alleles)[!allele_is_known], expected_unknown)
    J_aa <- translate_J_alleles(J_alleles, auxdata)
    expect_identical(unname(is.na(J_aa)), !allele_is_known)
    expect_identical(J_aa[["IGHJ1*01"]], "AEYFQHWGQGTLVTVSS")
    expect_true(all(grepl("[WF]G.G", J_aa[allele_is_known])))

    ## --- for mouse J alleles (from IMGT) ---

    auxdata <- load_auxdata("mouse", which="original")

    db_name <- "IMGT-202614-2.Mus_musculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    ## 3 mouse J alleles in IMGT release 202614-2 have no entries
    ## in 'auxdata'. Note that this could change in the future.
    expected_unknown <- c("IGLJ2P*01", "IGLJ4*01_Mus_spretus",
                          "IGLJ5*01_Mus_spretus")
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    expect_identical(names(J_alleles)[!allele_is_known], expected_unknown)
    ## Get rid of the "unknown" J alleles.
    known_J_alleles <- J_alleles[allele_is_known]
    J_aa <- translate_J_alleles(known_J_alleles, auxdata)
    expect_false(anyNA(J_aa))
    expect_identical(J_aa[["IGHJ1*01"]], "YWYFDVWGAGTTVTVSS")
})

test_that("J_allele_has_stop_codon()", {

    ## --- for human J alleles (from OGRDB and IMGT) ---

    auxdata <- load_and_fix_human_auxdata()

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    has_stop_codon <- J_allele_has_stop_codon(J_alleles, auxdata)
    expect_true(is.logical(has_stop_codon))
    expect_identical(names(has_stop_codon), names(J_alleles))
    expect_false(any(has_stop_codon))

    db_name <- "IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    has_stop_codon <- J_allele_has_stop_codon(J_alleles, auxdata)
    expect_identical(unname(is.na(has_stop_codon)), !allele_is_known)
    expect_false(any(has_stop_codon[allele_is_known]))

    ## --- for rabbit J alleles (from IMGT) ---

    auxdata <- load_auxdata("rabbit", which="original")

    db_name <- "IMGT-202614-2.Oryctolagus_cuniculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    has_stop_codon <- J_allele_has_stop_codon(J_alleles, auxdata)
    expect_false(anyNA(has_stop_codon))
    expect_identical(names(J_alleles)[has_stop_codon], "IGKJ1-2*04")
})

test_that("translate_fwr4()", {

    ## --- for human J alleles (from OGRDB and IMGT) ---

    auxdata <- load_and_fix_human_auxdata()

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    fwr4_aa <- translate_fwr4(J_alleles, auxdata)
    expect_true(is.character(fwr4_aa))
    expect_identical(names(fwr4_aa), names(J_alleles))
    expect_false(anyNA(fwr4_aa))
    fwr4_head <- translate_fwr4(J_alleles, auxdata, max.codons=4L)
    expect_true(is.character(fwr4_head))
    expect_identical(names(fwr4_head), names(J_alleles))
    expect_false(anyNA(fwr4_head))
    expect_true(all(nchar(fwr4_head) == 4L))
    expect_true(all(grepl("^[WF]G.G$", fwr4_head)))

    db_name <- "IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    fwr4_head <- translate_fwr4(J_alleles, auxdata, max.codons=4L)
    expect_identical(unname(is.na(fwr4_head)), !allele_is_known)
    expect_true(all(grepl("^[WF]G.G$", fwr4_head[allele_is_known])))

    ## --- for mouse J alleles (from IMGT) ---

    auxdata <- load_auxdata("mouse", which="original")

    db_name <- "IMGT-202614-2.Mus_musculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    ## 3 mouse J alleles in IMGT release 202614-2 have no entries
    ## in 'auxdata'. Note that this could change in the future.
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    expect_true(sum(!allele_is_known) <= 3L)
    ## Get rid of the "unknown" J alleles.
    known_J_alleles <- J_alleles[allele_is_known]
    fwr4_head <- translate_fwr4(known_J_alleles, auxdata, max.codons=4L)
    expect_false(anyNA(fwr4_head))
    ## 3 "known" mouse J alleles in IMGT release 202614-2 don't have
    ## the expected motif at the beginning of their FWR4 region.
    ## Is this expected? Could this change in the future?
    surprise <- fwr4_head[!grepl("^[WF]G.G$", fwr4_head)]
    expect_identical(names(surprise), c("IGKJ3*01", "IGKJ3*02", "IGLJ3P*01"))
    expect_identical(unname(surprise), c("FSDG", "FSDG", "FSSN"))

    ## --- for rat J alleles (from IMGT) ---

    auxdata <- load_auxdata("rat", which="original")

    db_name <- "IMGT-202614-2.Rattus_norvegicus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    fwr4_head <- translate_fwr4(J_alleles, auxdata, max.codons=4L)
    ## translate_fwr4() uses 'auxdata$cdr3_end' to get the position of
    ## the first FWR4 codons, but this column has an NA for IGKJ3*01.
    ## This could change in the future.
    ok <- !is.na(fwr4_head)
    expect_identical(names(fwr4_head)[!ok], "IGKJ3*01")
    ## Get rid of IGKJ3*01.
    fwr4_head <- fwr4_head[ok]
    ## 2 "known" rat J alleles in IMGT release 202614-2 don't have
    ## the expected motif at the beginning of their FWR4 region.
    ## Is this expected? Could this change in the future?
    surprise <- fwr4_head[!grepl("^[WF]G.G$", fwr4_head)]
    expect_identical(names(surprise), c("IGLJ2*01", "IGLJ4*01"))
    expect_identical(unname(surprise), c("LGKG", "LGKG"))
})


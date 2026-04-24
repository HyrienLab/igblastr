test_that("install_IMGT_germline_db()", {
    db_name <- install_IMGT_germline_db("202614-2", "Homo sapiens",
                                        overwrite=TRUE)
    expect_identical(db_name, "IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL")
    intdata <- load_intdata(db_name)
    igblastr:::check_ndm_data_col2class(intdata)
    V_allele_names <- names(load_germline_db(db_name, region_types="V"))
    expect_identical(intdata[ , "allele_name"], V_allele_names)
    auxdata <- load_auxdata(db_name)
    igblastr:::check_auxdata_col2class(auxdata)
    J_allele_names <- names(load_germline_db(db_name, region_types="J"))
    expect_identical(auxdata[ , "allele_name"], J_allele_names)
    use_germline_db(db_name)
    rm_germline_db(db_name)

    db_name <- install_IMGT_germline_db("202614-2", "Homo sapiens",
                                        tcr.db=TRUE, overwrite=TRUE)
    expect_identical(db_name, "IMGT-202614-2.Homo_sapiens.TRA+TRB+TRG+TRD")
    intdata <- load_intdata(db_name)
    igblastr:::check_ndm_data_col2class(intdata)
    V_allele_names <- names(load_germline_db(db_name, region_types="V"))
    expect_identical(intdata[ , "allele_name"], V_allele_names)
    use_germline_db(db_name)
    rm_germline_db(db_name)

    db_name <- install_IMGT_germline_db("202614-2", "Homo sapiens",
                                        loci="IGH", overwrite=TRUE)
    expect_identical(db_name, "IMGT-202614-2.Homo_sapiens.IGH")
    intdata <- load_intdata(db_name)
    igblastr:::check_ndm_data_col2class(intdata)
    V_allele_names <- names(load_germline_db(db_name, region_types="V"))
    expect_identical(intdata[ , "allele_name"], V_allele_names)
    auxdata <- load_auxdata(db_name)
    igblastr:::check_auxdata_col2class(auxdata)
    J_allele_names <- names(load_germline_db(db_name, region_types="J"))
    expect_identical(auxdata[ , "allele_name"], J_allele_names)
    use_germline_db(db_name)
    rm_germline_db(db_name)
 
    db_name <- install_IMGT_germline_db("202614-2", "Homo sapiens",
                                        loci=c("TRB", "TRA"), overwrite=TRUE)
    expect_identical(db_name, "IMGT-202614-2.Homo_sapiens.TRA+TRB")
    intdata <- load_intdata(db_name)
    igblastr:::check_ndm_data_col2class(intdata)
    V_allele_names <- names(load_germline_db(db_name, region_types="V"))
    expect_identical(intdata[ , "allele_name"], V_allele_names)
    use_germline_db(db_name)
    rm_germline_db(db_name)

    ## Yep, all the FASTA files for Mus_musculus_C57BL6J in IMGT release
    ## 202530-1 are empty!
    expect_error(
        install_IMGT_germline_db("202530-1", "Mus_musculus_C57BL6J",
                                 without.intdata=TRUE, without.auxdata=TRUE,
                                 overwrite=TRUE),
        regexp="no alleles found in FASTA files"
    )
    expect_error(
        install_IMGT_germline_db("202530-1", "Mus_musculus_C57BL6J",
                                 tcr.db=TRUE,
                                 without.intdata=TRUE, without.auxdata=TRUE,
                                 overwrite=TRUE),
        regexp="no alleles found in FASTA files"
    )
})

test_that("install_IMGT_c_region_db()", {
    errmsg <- "The following allele names are ambiguous:"
    expect_error2(
        install_IMGT_c_region_db("human", "IGH+IGK", overwrite=TRUE),
        errmsg
    )
    db_name <- install_IMGT_c_region_db("human", "IGH+IGK",
                                        disambiguate.allele.names=TRUE,
                                        overwrite=TRUE)
    expect_identical(db_name, "IMGT.human.IGH+IGK.202412")
    use_c_region_db(db_name)
    rm_c_region_db(db_name)

    db_name <- install_IMGT_c_region_db("Homo sapiens", "TRB+TRA",
                                        overwrite=TRUE)
    expect_identical(db_name, "IMGT.human.TRA+TRB.202509")
    use_c_region_db(db_name)
    rm_c_region_db(db_name)
})


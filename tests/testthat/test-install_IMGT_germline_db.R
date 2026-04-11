test_that("install_IMGT_germline_db()", {
    db_name <- install_IMGT_germline_db("202614-2", "Homo sapiens",
                                        overwrite=TRUE)
    expect_identical(db_name, "IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL")
    intdata <- load_intdata(db_name)
    igblastr:::check_ndm_data_col2class(intdata)
    use_germline_db(db_name)
    rm_germline_db(db_name)

    db_name <- install_IMGT_germline_db("202614-2", "Homo sapiens",
                                        tcr.db=TRUE, overwrite=TRUE)
    expect_identical(db_name, "IMGT-202614-2.Homo_sapiens.TRA+TRB+TRG+TRD")
    intdata <- load_intdata(db_name)
    igblastr:::check_ndm_data_col2class(intdata)
    use_germline_db(db_name)
    rm_germline_db(db_name)

    db_name <- install_IMGT_germline_db("202614-2", "Homo sapiens",
                                        loci="IGH", overwrite=TRUE)
    expect_identical(db_name, "IMGT-202614-2.Homo_sapiens.IGH")
    intdata <- load_intdata(db_name)
    igblastr:::check_ndm_data_col2class(intdata)
    use_germline_db(db_name)
    rm_germline_db(db_name)
 
    db_name <- install_IMGT_germline_db("202614-2", "Homo sapiens",
                                        loci=c("TRB", "TRA"), overwrite=TRUE)
    expect_identical(db_name, "IMGT-202614-2.Homo_sapiens.TRA+TRB")
    intdata <- load_intdata(db_name)
    igblastr:::check_ndm_data_col2class(intdata)
    use_germline_db(db_name)
    rm_germline_db(db_name)
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


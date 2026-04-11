### Install IMGT germline dbs used in tests below. Note that only the first
### installation actually triggers a download from IMGT. All subsequent
### installations obtain the data from the IMGT local store (located
### in 'igblastr_cache(IMGT_STORE)') so are very fast and work offline.
install_IMGT_germline_db("202614-2", "Oryctolagus cuniculus",
                         without.auxdata=TRUE, overwrite=TRUE)

.AUXDATA_COLNAMES <- names(igblastr:::AUXDATA_COL2CLASS)

test_that("load_auxdata()", {
    organisms <- list_igblast_organisms()
    for (organism in organisms) {
        auxdata <- load_auxdata(organism, which="original")
        expect_true(is.data.frame(auxdata))
        expect_identical(colnames(auxdata), .AUXDATA_COLNAMES)
        ## 1 row is repeated in human_gl.aux (the row for TRAJ13*02)
        if (organism == "human") {
            ok <- igblastr:::rows_with_same_key_are_identical(auxdata,
                                                              "allele_name")
            expect_true(ok)
        } else {
            expect_identical(anyDuplicated(auxdata[ , "allele_name"]), 0L)
        }
    }
})

test_that("compute_germline_db_auxdata()", {
    ## Compare computed auxdata for
    ## IMGT-202614-2.Oryctolagus_cuniculus.IGH+IGK+IGL
    ## with rabbit auxdata included in IgBLAST.

    db_name <- "IMGT-202614-2.Oryctolagus_cuniculus.IGH+IGK+IGL"
    auxdata <- compute_germline_db_auxdata(db_name)
    auxdata0 <- load_auxdata("rabbit")
    m <- match(auxdata[ , "allele_name"], auxdata0[ , "allele_name"])
    ## All J alleles in germline db are annotated in IgBLAST auxdata.
    expect_false(anyNA(m))
    ## Check that NA-free rows in 'auxdata' contain the same information as
    ## their corresponding rows in 'auxdata0'.
    auxdata1 <- S4Vectors:::extract_data_frame_rows(auxdata0, m)
    keep_idx <- which(!is.na(auxdata[ , "cdr3_end"]))
    expect_identical(auxdata[keep_idx, ], auxdata1[keep_idx, ])
})


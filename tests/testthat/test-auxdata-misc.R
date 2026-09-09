### Install IMGT germline dbs used in tests below. Note that only the first
### installation actually triggers a download from IMGT. All subsequent
### installations obtain the data from the IMGT local store (located
### in 'igblastr_cache(IMGT_STORE)') so are very fast and work offline.
install_IMGT_germline_db("202631-1", "Oryctolagus cuniculus",
                         auto.auxdata=FALSE, overwrite=TRUE)

test_that("compute_germline_db_auxdata()", {
    ## Compare computed auxdata for
    ## IMGT-202631-1.Oryctolagus_cuniculus.IGH+IGK+IGL
    ## with rabbit auxdata included in IgBLAST.

    db_name <- "IMGT-202631-1.Oryctolagus_cuniculus.IGH+IGK+IGL"
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


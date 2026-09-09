test_that("list_germline_dbs()", {
    df <- list_germline_dbs()
    expect_true(is.data.frame(df))
    expected_colnames <- c("db_name", "V", "D", "J", "intdata", "auxdata")
    expect_identical(colnames(df), expected_colnames)

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
    use_germline_db(db_name)
    printed <- print(df)
    expect_true(is.data.frame(printed))
    expect_equal(nrow(printed), nrow(df))
    ## One extra column for the asterisk.
    expect_equal(ncol(printed), ncol(df) + 1L)
    expect_identical(trimws(colnames(printed)), c(expected_colnames, ""))

    ## Check consistency of counts reported by short and long listings.
    install_IMGT_germline_db("202631-1", "Homo sapiens", overwrite=TRUE)
    df <- list_germline_dbs()  # short listing
    all_counts <- list_germline_dbs(long.listing=TRUE)  # long listing
    expect_identical(nrow(df), length(all_counts))
    for (i in seq_along(all_counts)) {
        counts <- all_counts[[i]]
        for (region_type in c("V", "D", "J"))
            expect_identical(sum(counts[ , region_type]), df[i , region_type])
    }
})


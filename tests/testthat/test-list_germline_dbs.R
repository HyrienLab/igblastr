test_that("use_germline_db()", {
    db_name <- "_AIRR.human.IGH+IGK+IGL.202501"
    use_germline_db(db_name)
    expect_identical(use_germline_db(), db_name)
})

test_that("list_germline_dbs()", {
    df <- list_germline_dbs()
    expect_true(is.data.frame(df))
    expected_colnames <- c("db_name",
                           "IGHV", "IGHD", "IGHJ",
                           "IGKV", "IGKJ", "IGLV", "IGLJ")
    expect_identical(colnames(df), expected_colnames)

    db_name <- "_AIRR.human.IGH+IGK+IGL.202501"
    use_germline_db(db_name)
    printed <- print(df)
    expect_true(is.data.frame(printed))
    expect_equal(nrow(printed), nrow(df))
    ## One extra column for the asterisk.
    expect_equal(ncol(printed), ncol(df) + 1L)
    expect_identical(trimws(colnames(printed)), c(expected_colnames, ""))
})


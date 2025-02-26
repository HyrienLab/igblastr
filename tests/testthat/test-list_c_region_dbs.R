test_that("use_c_region_db()", {
    use_c_region_db("")
    expect_identical(use_c_region_db(), "")

    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    expect_identical(use_c_region_db(), db_name)
})

test_that("list_c_region_dbs()", {
    df <- list_c_region_dbs()
    expect_true(is.data.frame(df))
    expected_colnames <- c("db_name", "IGH", "IGK", "IGL")
    expect_identical(colnames(df), expected_colnames)

    use_c_region_db("")
    printed <- print(df)
    expect_true(is.data.frame(printed))
    expect_identical(dim(printed), dim(df))

    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    printed <- print(df)
    expect_true(is.data.frame(printed))
    expect_equal(nrow(printed), nrow(df))
    ## One extra column for the asterisk.
    expect_equal(ncol(printed), ncol(df) + 1L)
    expect_identical(trimws(colnames(printed)), c(expected_colnames, ""))
})


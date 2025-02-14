test_that("use_germline_db()", {
    db_name <- "_AIRR.human.IGH+IGK+IGL.202501"
    use_germline_db(db_name)
    expect_identical(use_germline_db(), db_name)
})

test_that("list_germline_dbs()", {
    object <- list_germline_dbs()
    expect_true(is.data.frame(object))
    expected_colnames <- c("db_name",
                           "IGHV", "IGHD", "IGHJ",
                           "IGKV", "IGKJ", "IGLV", "IGLJ")
    expect_identical(colnames(object), expected_colnames)

    db_name <- "_AIRR.human.IGH+IGK+IGL.202501"
    use_germline_db(db_name)
    printed <- print(object)
    expect_true(is.data.frame(printed))
    expect_equal(nrow(printed), nrow(object))
    ## One extra column for the asterisk.
    expect_equal(ncol(printed), ncol(object) + 1L)
    expect_identical(trimws(colnames(printed)), c(expected_colnames, ""))
})


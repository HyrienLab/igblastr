test_that("list_germline_dbs()", {
    object <- list_germline_dbs()
    expect_true(is.data.frame(object))
    expected_colnames <- c("db_name",
                           "IGHV", "IGHD", "IGHJ",
                           "IGKV", "IGKJ", "IGLV", "IGLJ")
    expect_identical(colnames(object), expected_colnames)
})

test_that("use_germline_db()", {
    db_name <- "_AIRR.human.IGH+IGK+IGL.202412"
    use_germline_db(db_name)
    expect_identical(use_germline_db(), db_name)
    object <- list_germline_dbs()
    expect_true(is.data.frame(object))
    expected_colnames <- c("db_name",
                           "IGHV", "IGHD", "IGHJ",
                           "IGKV", "IGKJ", "IGLV", "IGLJ", " ")
    expect_identical(colnames(object), expected_colnames)
})


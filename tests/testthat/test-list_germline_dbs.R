test_that("list_germline_dbs()", {
    object <- list_germline_dbs()
    expect_true(is.data.frame(object))
    expected_colnames <- c("db_name",
                           "IGHV", "IGHD", "IGHJ",
                           "IGKV", "IGKJ", "IGLV", "IGLJ",
                           "used")
    expect_identical(colnames(object), expected_colnames)
})


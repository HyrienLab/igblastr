test_that("list_germline_dbs()", {
    object <- list_germline_dbs()
    expect_true(is.data.frame(object))
    expected_colnames <- c("db_name",
                           "IGHV", "IGHD", "IGHJ",
                           "IGKV", "IGKJ", "IGLV", "IGLJ",
                           " ")
    expect_identical(colnames(object), expected_colnames)
})


test_that("list_c_region_dbs()", {
    object <- list_c_region_dbs()
    expect_true(is.data.frame(object))
    expect_identical(colnames(object), c("db_name", "nCregions", "used"))
})

test_that("use_c_region_db()", {
    use_c_region_db("")
    expect_identical(use_c_region_db(), "")
    use_c_region_db("IMGT.rabbit.IG")
    expect_identical(use_c_region_db(), "IMGT.rabbit.IG")
})


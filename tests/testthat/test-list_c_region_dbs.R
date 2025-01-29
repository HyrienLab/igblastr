test_that("list_c_region_dbs()", {
    use_c_region_db("")
    object <- list_c_region_dbs()
    expect_true(is.data.frame(object))
    expected <- c("db_name", "IGH", "IGK", "IGL")
    expect_identical(colnames(object), expected)
})

test_that("use_c_region_db()", {
    use_c_region_db("")
    expect_identical(use_c_region_db(), "")
    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    expect_identical(use_c_region_db(), db_name)
    object <- list_c_region_dbs()
    expect_true(is.data.frame(object))
    expected <- c("db_name", "IGH", "IGK", "IGL", " ")
    expect_identical(colnames(object), expected)
})


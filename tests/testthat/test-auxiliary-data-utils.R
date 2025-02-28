test_that("load_igblast_auxiliary_data()", {
    organisms <- list_igblast_organisms()
    for (organism in organisms) {
        df <- load_igblast_auxiliary_data(organism)
        expect_true(is.data.frame(df))
        expect_equal(ncol(df), 5)
    }
})


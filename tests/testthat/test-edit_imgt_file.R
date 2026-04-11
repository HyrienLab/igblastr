test_that("redit_imgt_file() behaves like edit_imgt_file()", {
    ## Requires Perl!
    failures <- igblastr:::validate_redit_imgt_file_on_IMGT_release("202614-2")
    expect_identical(failures, 0L)
})


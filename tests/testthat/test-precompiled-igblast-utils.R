test_that("inst/extdata/ncbi_igblast_data_files/ is not out-of-sync", {
    ## igblastr:::check_ncbi_igblast_data_files() is not guaranteed to work
    ## on Windows.
    OS_arch <- igblastr:::get_OS_arch()
    OS <- tolower(OS_arch[["OS"]])
    if (OS != "windows") {
	prev_warn <- getOption("warn")
        options(warn=2)
        on.exit(options(warn=prev_warn))
        igblastr:::check_ncbi_igblast_data_files(get_igblast_root())
    }
})


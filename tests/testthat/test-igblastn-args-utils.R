test_that("make_igblastn_command_line_args()", {
    ## make_igblastn_command_line_args() is not exported.
    make_igblastn_command_line_args <-
        igblastr:::make_igblastn_command_line_args
    CORE_COLNAMES <- c("query", "outfmt", "organism",
                       paste0("germline_db_", c("V", "D", "J")))

    ## For this test, we use invalid V-, D-, J-region dbs but
    ## make_igblastn_command_line_args() should still accept them.
    ## Should always work, even if no germline db is currently in use.
    db_path <- system.file(package="igblastr", "extdata")
    germline_db_V <- file.path(db_path, "V")
    germline_db_D <- file.path(db_path, "D")
    germline_db_J <- file.path(db_path, "J")
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                     organism="rhesus_monkey",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     c_region_db=NULL)
    expect_true(is.character(cmd_args))
    expect_identical(names(cmd_args), CORE_COLNAMES)

    ## For this test we're not specifying any of the 'germline_db_[VDJ]'
    ## argument so we need to select a local germline db. Once we do
    ## this, we don't need to specify 'organism' either.
    db_name <- "_AIRR.human.IGH+IGK+IGL.202501"
    use_germline_db(db_name)
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                c_region_db=NULL)
    expect_true(is.character(cmd_args))
    expect_identical(names(cmd_args), CORE_COLNAMES)

    use_c_region_db("")
    cmd_args <- make_igblastn_command_line_args("path/to/query")
    expect_true(is.character(cmd_args))
    expect_identical(names(cmd_args), CORE_COLNAMES)

    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    cmd_args <- make_igblastn_command_line_args("path/to/query")
    expected_colnames <- c(CORE_COLNAMES, "c_region_db")
    expect_identical(names(cmd_args), expected_colnames)
})


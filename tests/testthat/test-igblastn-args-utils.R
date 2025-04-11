test_that("make_igblastn_command_line_args()", {
    ## make_igblastn_command_line_args() is not exported.
    make_igblastn_command_line_args <-
        igblastr:::make_igblastn_command_line_args
    CORE_ARGNAMES <- c("query", "outfmt",
                       paste0("germline_db_", c("V", "D", "J")),
                       "organism")

    ## For this test, we use invalid V-, D-, J-region dbs but
    ## make_igblastn_command_line_args() should still accept them.
    ## Should always work, even if no germline db is currently in use.
    db_path <- system.file(package="igblastr", "extdata")
    germline_db_V <- file.path(db_path, "V")
    germline_db_D <- file.path(db_path, "D")
    germline_db_J <- file.path(db_path, "J")
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     organism="rhesus_monkey",
                     c_region_db=NULL,
                     auxiliary_data=NULL)
    expect_true(is.character(cmd_args))
    expect_identical(names(cmd_args), CORE_ARGNAMES)

    ## For this test we're not specifying any of the 'germline_db_[VDJ]'
    ## arguments so we need to select a local germline db. Once we do
    ## this, we don't need to specify 'organism' either.
    db_name <- "_AIRR.human.IGH+IGK+IGL.202501"
    use_germline_db(db_name)
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                c_region_db=NULL,
                                                auxiliary_data=NULL)
    expect_true(is.character(cmd_args))
    expect_identical(names(cmd_args), CORE_ARGNAMES)

    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                c_region_db=NULL)
    expect_true(is.character(cmd_args))
    expected_argnames <- c(CORE_ARGNAMES, "auxiliary_data")
    expect_identical(names(cmd_args), expected_argnames)

    use_c_region_db("")
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                auxiliary_data=NULL)
    expect_true(is.character(cmd_args))
    expect_identical(names(cmd_args), CORE_ARGNAMES)

    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                auxiliary_data=NULL)
    expect_true(is.character(cmd_args))
    expected_argnames <- c(CORE_ARGNAMES, "c_region_db")
    expect_identical(names(cmd_args), expected_argnames)

    cmd_args <- make_igblastn_command_line_args("path/to/query")
    expect_true(is.character(cmd_args))
    expected_argnames <- c(CORE_ARGNAMES, "c_region_db", "auxiliary_data")
    expect_identical(names(cmd_args), expected_argnames)
})


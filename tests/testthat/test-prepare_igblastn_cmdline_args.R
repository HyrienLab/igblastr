test_that(".normarg_custom_internal_data()", {
    normarg_custom_internal_data <- igblastr:::.normarg_custom_internal_data

    ## Note that we only check "auto" mode for now.

    ## --- With a supplied 'germline_db_V' that looks like the path ---
    ## --- to a V "blast db" that belongs to a cached germline db.  ---

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"  # includes internal data
    db_path <- igblastr:::get_germline_db_path(db_name)
    germline_db_V <- file.path(db_path, "V")

    custom_internal_data <- normarg_custom_internal_data("auto",
                                    germline_db_V, domain_system="imgt")
    expected <- file.path(db_path, "internal_data", "V.ndm.imgt")
    expect_identical(custom_internal_data, expected)
    expect_true(file.exists(custom_internal_data))

    ## Returns NULL with a warning.
    regexp <- paste0("internal data file .*\\.kabat not found.*",
                     "--> using IgBLAST internal data")
    expect_warning(
        custom_internal_data <- normarg_custom_internal_data("auto",
                                        germline_db_V, domain_system="kabat"),
        regexp=regexp
    )
    expect_true(is.null(custom_internal_data))

    ## No internal data.
    db_name <- install_IMGT_germline_db("202518-3", "Homo sapiens",
                                        without.intdata=TRUE, overwrite=TRUE)
    db_path <- igblastr:::get_germline_db_path(db_name)
    germline_db_V <- file.path(db_path, "V")

    ## Returns NULL with no warning.
    custom_internal_data <- normarg_custom_internal_data("auto",
                                    germline_db_V, domain_system="imgt")
    expect_true(is.null(custom_internal_data))
    custom_internal_data <- normarg_custom_internal_data("auto",
                                    germline_db_V, domain_system="kabat")
    expect_true(is.null(custom_internal_data))

    ## --- With a supplied 'germline_db_V' that does NOT look like the  ---
    ## --- path to a V "blast db" that belongs to a cached germline db. ---

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"  # includes internal data
    db_path <- igblastr:::get_germline_db_path(db_name)
    germline_db_V <- file.path(db_path, "J")
    custom_internal_data <- normarg_custom_internal_data("auto",
                                    germline_db_V, domain_system="imgt")
    expect_true(is.null(custom_internal_data))

    germline_db_V <- file.path(tempfile(), "V")
    custom_internal_data <- normarg_custom_internal_data("auto",
                                    germline_db_V, domain_system="imgt")
    expect_true(is.null(custom_internal_data))
})

test_that(".normarg_auxiliary_data()", {
    normarg_auxiliary_data <- igblastr:::.normarg_auxiliary_data

    expect_error(normarg_auxiliary_data(2), regexp="must be")
    expect_true(is.null(normarg_auxiliary_data(NULL)))
    expect_error(normarg_auxiliary_data(data.frame(aa=1:3)))

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    use_germline_db(db_name)
    db_path <- igblastr:::get_germline_db_path(db_name)
    germline_db_J <- file.path(db_path, "J")

    auxdata_path <- get_auxdata_path("human")
    auxiliary_data <- normarg_auxiliary_data(auxdata_path,
                                             germline_db_J=germline_db_J)
    expect_identical(auxiliary_data, auxdata_path)

    auxdata <- load_auxdata("human")
    auxiliary_data <- normarg_auxiliary_data(auxdata,
                                             germline_db_J=germline_db_J)
    expect_identical(read_auxdata(auxiliary_data), auxdata)
    expect_identical(attr(auxiliary_data, "safe_to_remove"), TRUE)

    expect_error(
        normarg_auxiliary_data(rev(auxdata)),
        regexp="must have the following columns"
    )

    auxdata_path <- get_auxdata_path("mouse")
    expect_warning(
        auxiliary_data <- normarg_auxiliary_data(auxdata_path,
                                                 germline_db_J=germline_db_J),
        regexp="Incomplete auxiliary data"
    )
    expect_identical(auxiliary_data, auxdata_path)

    auxdata <- load_auxdata("mouse")
    expect_warning(
        auxiliary_data <- normarg_auxiliary_data(auxdata,
                                                 germline_db_J=germline_db_J),
        regexp="Incomplete auxiliary data"
    )
    expect_identical(read_auxdata(auxiliary_data), auxdata)
    expect_identical(attr(auxiliary_data, "safe_to_remove"), TRUE)

    auxiliary_data <- normarg_auxiliary_data("auto", organism="human",
                                             germline_db_J=germline_db_J)
    expect_identical(auxiliary_data, get_auxdata_path(db_name))

    auxiliary_data <- normarg_auxiliary_data("auto", organism=NULL,
                                             germline_db_J=germline_db_J,
                                             no_auto_germline_dbs=FALSE)
    expect_identical(auxiliary_data, get_auxdata_path(db_name))

    db_name <- "_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202501"
    use_germline_db(db_name)
    db_path <- igblastr:::get_germline_db_path(db_name)
    germline_db_J <- file.path(db_path, "J")

    auxiliary_data <- normarg_auxiliary_data("auto", organism=NULL,
                                             germline_db_J=germline_db_J,
                                             no_auto_germline_dbs=FALSE)
    expect_identical(auxiliary_data, get_auxdata_path("mouse"))

    expect_error(
        normarg_auxiliary_data("auto", organism=NULL,
                               germline_db_J=germline_db_J,
                               no_auto_germline_dbs=TRUE),
        regexp="Failed to automatically figure out"
    )
})

test_that("prepare_igblastn_cmdline_args()", {
    ## prepare_igblastn_cmdline_args() is not exported.
    prepare_igblastn_cmdline_args <- igblastr:::prepare_igblastn_cmdline_args
    CORE_ARGNAMES <- c("query", "outfmt",
                       paste0("germline_db_", c("V", "D", "J")),
                       "organism", "domain_system", "ig_seqtype", "out")
    organism_idx <- match("organism", CORE_ARGNAMES)
    use_c_region_db("")

    ## --- no germline db selection ---

    ## For this test, we use invalid V-, D-, J-region dbs but
    ## prepare_igblastn_cmdline_args() should still accept them.
    ## Should always work, even if no germline db is currently in use.
    db_path <- system.file(package="igblastr", "extdata")
    germline_db_V <- file.path(db_path, "V")
    germline_db_D <- file.path(db_path, "D")
    germline_db_J <- file.path(db_path, "J")
    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     organism="rhesus_monkey",
                     c_region_db=NULL,
                     custom_internal_data=NULL,
                     auxiliary_data=NULL,
                     ig_seqtype="Ig")
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expect_identical(names(cmd_args), CORE_ARGNAMES)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$organism, "rhesus_monkey")
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    ## Error: 'organism' must be specified
    expect_error(prepare_igblastn_cmdline_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J),
                 regexp="organism.*must be specified")
    ## Error: 'ig_seqtype' must be specified
    expect_error(prepare_igblastn_cmdline_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     organism="rabbit"),
                 regexp="ig_seqtype.*must be specified")

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     organism="rabbit",
                     ig_seqtype="tcr")
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    opt_argnames <- "auxiliary_data"
    expected_argnames <- append(CORE_ARGNAMES, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    ## --- selecting _OGRDB.human.IGH+IGK+IGL.202410 ---

    ## For these tests we're not specifying any of the 'germline_db_[VDJ]'
    ## arguments so we need to select a cached germline db. Once we do
    ## this, we don't need to specify 'organism' either.

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    suppressMessages(use_germline_db(db_name))
    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              c_region_db=NULL,
                                              custom_internal_data=NULL,
                                              auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expect_identical(names(cmd_args), CORE_ARGNAMES)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$organism, "human")
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              c_region_db=NULL,
                                              auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$custom_internal_data, get_intdata_path(db_name))
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              c_region_db=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    opt_argnames <- "auxiliary_data"
    expected_argnames <- append(expected_argnames, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$custom_internal_data, get_intdata_path(db_name))
    expect_identical(cmd_args$auxiliary_data, get_auxdata_path(db_name))
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)

    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "c_region_db"
    opt_argnames <- "custom_internal_data"
    expected_argnames <- append(expected_argnames, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query")
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "c_region_db"
    opt_argnames <- c("custom_internal_data", "auxiliary_data")
    expected_argnames <- append(expected_argnames, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    ## --- selecting IMGT-202531-1.Mus_musculus.IGH+IGK+IGL ---

    db_name <- install_IMGT_germline_db("202531-1", "Mus musculus",
                                        tcr.db=TRUE, overwrite=TRUE)
    suppressMessages(use_germline_db(db_name))

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              c_region_db=NULL,
                                              auxiliary_data=NULL)
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "TCR")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    ## --- selecting IMGT-202518-3.Sus_scrofa.IGH+IGK+IGL ---
    ## Note that pig is not an IgBLAST organism!

    db_name <- install_IMGT_germline_db("202518-3", "Sus_scrofa",
                                        overwrite=TRUE)
    suppressMessages(use_germline_db(db_name))
    use_c_region_db("")

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              auxiliary_data=NULL)
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)

    db_name <- install_IMGT_germline_db("202518-3", "Sus_scrofa",
                                        without.intdata=TRUE,
                                        without.auxdata=TRUE, overwrite=TRUE)
    suppressMessages(use_germline_db(db_name))

    errmsg <- "Don't know how to infer 'organism' from germline db name"
    expect_error2(prepare_igblastn_cmdline_args("path/to/query"), errmsg)

    custom_internal_data <- get_intdata_path("_OGRDB.human.IGH+IGK+IGL.202410")
    regexp <- "what auxiliary data to use"
    expect_error(
        suppressWarnings(
            prepare_igblastn_cmdline_args(
                    "path/to/query",
                    custom_internal_data=custom_internal_data)
        ),
        regexp=regexp
    )

    expect_warning(
        cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                 custom_internal_data=custom_internal_data,
                                 auxiliary_data=NULL),
        regexp="Incomplete custom internal data"
    )
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)
})


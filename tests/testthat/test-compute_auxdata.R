### Install IMGT germline dbs used in tests below. Note that only the first
### installation actually triggers a download from IMGT. All subsequent
### installations obtain the data from the IMGT local store (located
### in 'igblastr_cache(IMGT_STORE)') so are very fast and work offline.
install_IMGT_germline_db("202614-2", "Homo sapiens",
                         without.auxdata=TRUE, overwrite=TRUE)
install_IMGT_germline_db("202614-2", "Mus musculus",
                         without.auxdata=TRUE, overwrite=TRUE)
install_IMGT_germline_db("202614-2", "Rattus norvegicus",
                         without.auxdata=TRUE, overwrite=TRUE)

.AUXDATA_COLNAMES <- names(igblastr:::AUXDATA_COL2CLASS)

test_that("compute_auxdata()", {

    ## --- for human J alleles (from AIRR and IMGT) ---

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- compute_auxdata(J_alleles)
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Now we're going to check that 'computed_auxdata' agrees with the
    ## auxiliary data included in IgBLAST. More precisely, we're going to
    ## check that it's a subset of 'load_auxdata("human", which="original")'.

    ## All the J alleles in _OGRDB.human.IGH+IGK+IGL.202410 are annotated
    ## in human_gl.aux so we expect no NAs in 'm' below.
    orig_auxdata <- load_and_fix_human_auxdata()
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    expect_false(anyNA(m))
    orig_auxdata <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m)
    expect_identical(computed_auxdata, orig_auxdata)

    db_name <- "IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- compute_auxdata(J_alleles)
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL
    ## are annotated in human_gl.aux so we expect a few NAs in 'm' below.
    orig_auxdata <- load_and_fix_human_auxdata()
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!is.na(m))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)

    ## --- for mouse J alleles (from IMGT) ---

    db_name <- "IMGT-202614-2.Mus_musculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- suppressWarnings(compute_auxdata(J_alleles))
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202614-2.Mus_musculus.IGH+IGK+IGL
    ## are annotated in mouse_gl.aux so we expect a few NAs in 'm' below.
    ## We will also skip validation for alleles for which no CDR3 end was
    ## found.
    orig_auxdata <- load_auxdata("mouse", which="original")
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!(is.na(computed_auxdata[ , "cdr3_end"]) | is.na(m)))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)

    ## --- for rat J alleles (from IMGT) ---

    db_name <- "IMGT-202614-2.Rattus_norvegicus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- suppressWarnings(compute_auxdata(J_alleles))
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202614-2.Mus_musculus.IGH+IGK+IGL
    ## are annotated in rat_gl.aux so we expect a few NAs in 'm' below.
    ## We will also skip validation for alleles for which no CDR3 end was
    ## found.
    orig_auxdata <- load_auxdata("rat", which="original")
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!(is.na(computed_auxdata[ , "cdr3_end"]) | is.na(m)))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)
})


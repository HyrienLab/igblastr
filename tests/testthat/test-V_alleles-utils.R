
test_that("translate_V_alleles()", {
    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    intdata <- load_intdata(db_name)
    V_alleles <- load_germline_db(db_name, region_types="V")

    V_aa <- translate_V_alleles(V_alleles, intdata)
    expect_true(is.character(V_aa))
    expect_identical(names(V_aa), names(V_alleles))
    expect_false(anyNA(V_aa))

    fwr1_aa <- translate_V_alleles(V_alleles, intdata, V_segment="fwr1")
    cdr1_aa <- translate_V_alleles(V_alleles, intdata, V_segment="cdr1")
    fwr2_aa <- translate_V_alleles(V_alleles, intdata, V_segment="fwr2")
    cdr2_aa <- translate_V_alleles(V_alleles, intdata, V_segment="cdr2")
    fwr3_aa <- translate_V_alleles(V_alleles, intdata, V_segment="fwr3")
    expect_error(translate_V_alleles(V_alleles, intdata, V_segment="cdr3"),
                 regexp="must be one of")

    all_fwrcdr_aa <- paste0(fwr1_aa, cdr1_aa, fwr2_aa, cdr2_aa, fwr3_aa)
    expect_identical(all_fwrcdr_aa,
                     unname(substr(V_aa, 1L, nchar(all_fwrcdr_aa))))
})

test_that("V_allele_has_stop_codon()", {
    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    intdata <- load_intdata(db_name)
    V_alleles <- load_germline_db(db_name, region_types="V")

    has_stop_codon <- V_allele_has_stop_codon(V_alleles, intdata)
    expect_true(is.logical(has_stop_codon))
    expect_identical(names(has_stop_codon), names(V_alleles))
    expect_false(anyNA(has_stop_codon))
})


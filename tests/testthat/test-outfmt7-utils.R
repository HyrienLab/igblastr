.check_fmt7record <- function(rec, expected_qseqid)
{
    expect_equal(class(rec), "fmt7record")
    expect_true(is.list(rec))
    expect_true(length(rec) %in% 5:6)
    expected_names <- c("query_details", "VDJ_rearrangement_summary",
                        "VDJ_junction_details",
                        "alignment_summary", "hit_table")
    if (length(rec) == 6L)
        expected_names <- append(expected_names, "subregion_sequence_details",
                                 after=3L)
    expect_identical(names(rec), expected_names)
    expect_equal(qseqid(rec), expected_qseqid)

    expect_identical(class(rec$query_details), "query_details")
    expect_identical(class(rec$VDJ_rearrangement_summary),
                     "VDJ_rearrangement_summary")
    expect_identical(class(rec$VDJ_junction_details), "VDJ_junction_details")
    if (length(rec) == 6L)
        expect_identical(class(rec$subregion_sequence_details),
                         "subregion_sequence_details")
    expect_identical(class(rec$alignment_summary),
                     c("alignment_summary", "data.frame"))
    expect_identical(class(rec$hit_table), c("hit_table", "data.frame"))
}

test_that("parse_outfmt7()", {
    catnap_bnabs <- system.file(package="igblastr", "extdata",
                                "catnap_bnabs.fasta")
    query <- head(readDNAStringSet(catnap_bnabs), n=10)
    db_name <- "_AIRR.human.IGH+IGK+IGL.202501"
    use_germline_db(db_name)

    ## parse_outfmt7() gets called internally by igblastn(..., outfmt=7).

    ## With 10 sequences.
    out <- igblastn(query, outfmt=7)
    expect_true(is.list(out))
    expect_equal(length(out), 2)
    expect_identical(names(out), c("records", "footer"))
    records <- out$records
    footer <- out$footer
    expect_true(is.list(records))
    expect_equal(length(records), length(query))
    expect_equal(class(footer), "fmt7footer")
    expect_true(is.character(footer))
    for (i in seq_along(records))
        .check_fmt7record(records[[i]], names(query)[[i]])

    ## With zero query sequences.
    out <- igblastn(query[0], outfmt=7)
    expect_true(is.list(out))
    expect_equal(length(out), 2)
    expect_identical(names(out), c("records", "footer"))
    records <- out$records
    footer <- out$footer
    expect_true(is.list(records))
    expect_equal(length(records), 0)
    expect_equal(class(footer), "fmt7footer")
    expect_true(is.character(footer))
})


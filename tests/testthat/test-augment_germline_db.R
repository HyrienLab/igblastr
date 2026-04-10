### The augment_germline_db_*() functions need to be revisited so
### that they do the right thing with internal and auxiliary data
### so we disable this for now.
if (FALSE) {
test_that("augment_germline_db_*()", {
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    query <- readDNAStringSet(query)[c(1:2, 117:120)]

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    use_germline_db(db_name)
    AIRR_df0 <- igblastn(query)

    destdir <- tempfile()

    ## Add 0 novel J gene alleles to 'db_name'.

    novel_J_alleles <- DNAStringSet()
    augment_germline_db_J(db_name, novel_J_alleles, destdir=destdir)
    all_J_alleles <- readDNAStringSet(file.path(destdir, "J.fasta"))
    expected <- load_germline_db(db_name, "J")
    expect_identical(as.character(all_J_alleles), as.character(expected))

    AIRR_df1 <- igblastn(query, germline_db_J=destdir)
    expect_identical(AIRR_df0, AIRR_df1)
    expect_error(augment_germline_db_J(db_name, novel_J_alleles,
                                       destdir=destdir))
    augment_germline_db_J(db_name, novel_J_alleles, destdir=destdir,
                          overwrite=TRUE)
    expect_error(igblastn(query, germline_db_V=destdir))

    ## Add 3 made-up novel V gene alleles to 'db_name'.

    novel_V_alleles <- system.file(package="igblastr", "extdata",
                                   "novel_germline_alleles",
                                   "fake_human_V_alleles.fasta")
    augment_germline_db_V(db_name, novel_V_alleles, destdir=destdir)
    AIRR_df2 <- igblastn(query, germline_db_V=destdir)
    idx <- which(AIRR_df0$v_call != AIRR_df2$v_call)
    expect_identical(idx, 2:3)
    expect_true(all.equal(AIRR_df0[-idx, ], AIRR_df2[-idx, ]))

    AIRR_df3 <- igblastn(query, germline_db_V=destdir, germline_db_J=destdir)
    expect_identical(AIRR_df2, AIRR_df3)
})
}

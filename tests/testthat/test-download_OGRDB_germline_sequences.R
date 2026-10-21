.TEST_GERMLINE_SETS <- list(

    list(organism="Homo sapiens",
         germline_sets=c(IGH_VDJ= 7, IGKappa_VJ=2, IGLambda_VJ=1),
         destfiles_per_set=c(3L, 2L, 2L)),

    list(organism="Homo sapiens",
         germline_sets=c(IGH_VDJ= 8, IGKappa_VJ=3, IGLambda_VJ=2),
         destfiles_per_set=c(3L, 2L, 2L)),

    list(organism="Homo sapiens",
         germline_sets=c(IGH_VDJ= 9, IGKappa_VJ=4, IGLambda_VJ=3),
         destfiles_per_set=c(3L, 2L, 2L)),

    list(organism="Homo sapiens",
         germline_sets=c(IGH_VDJ=10, IGKappa_VJ=5, IGLambda_VJ=4),
         destfiles_per_set=c(3L, 2L, 2L)),

    list(organism="Mus musculus",
         germline_sets=c(`C57BL/6 IGH`=5,
                         `C57BL/6J IGKV`=1, `IGKJ (all strains)`=1,
                         `C57BL/6J IGLV`=1, `IGLJ (all strains)`=1),
         destfiles_per_set=c(3L, 1L, 1L, 1L, 1L)),

    list(organism="Mus musculus",
         germline_sets=c(`CAST/EiJ IGH`=2,
                         `CAST/EiJ IGKV`=1, `IGKJ (all strains)`=1,
                         `CAST/EiJ IGLV`=2, `IGLJ (all strains)`=1),
         destfiles_per_set=c(3L, 1L, 1L, 1L, 1L)),

    list(organism="Mus musculus",
         germline_sets=c(`LEWES/EiJ IGH`=2,
                         `LEWES/EiJ IGKV`=1, `IGKJ (all strains)`=1,
                         `LEWES/EiJ IGLV`=1, `IGLJ (all strains)`=1),
         destfiles_per_set=c(3L, 1L, 1L, 1L, 1L)),

    list(organism="Mus musculus",
         germline_sets=c(`MSM/MsJ IGH`=1,
                         `MSM/MsJ IGKV`=1, `IGKJ (all strains)`=1,
                         `MSM/MsJ IGLV`=2, `IGLJ (all strains)`=1),
         destfiles_per_set=c(3L, 1L, 1L, 1L, 1L)),

    list(organism="Mus musculus",
         germline_sets=c(`NOD/ShiLtJ IGH`=1,
                         `NOD/ShiLtJ IGKV`=1, `IGKJ (all strains)`=1,
                         `NOD/ShiLtJ IGLV`=1, `IGLJ (all strains)`=1),
         destfiles_per_set=c(3L, 1L, 1L, 1L, 1L)),

    list(organism="Mus musculus",
         germline_sets=c(`PWD/PhJ IGH`=2,
                         `PWD/PhJ IGKV`=1, `IGKJ (all strains)`=1,
                         `PWD/PhJ IGLV`=1, `IGLJ (all strains)`=1),
         destfiles_per_set=c(3L, 1L, 1L, 1L, 1L)),

    list(organism="Macaca mulatta",
         germline_sets=c(IGH_VDJ=2, IGK_VJ=2, IGL_VJ=2),
         destfiles_per_set=c(3L, 2L, 2L)),

    ## This germline set does not pass the sanity checks performed by
    ## extract_alleles_from_ogrdb_json().
    ## Reported to William Lees by email on Oct 21, 2026.
    #list(organism="Oncorhynchus mykiss",
    #     germline_sets=c(IGH_VDJ=1),
    #     destfiles_per_set=3L),

    list(organism="Peromyscus leucopus",
         germline_sets=c(IGH_VDJ=2, IGK_VJ=2, IGL_VJ=1),
         destfiles_per_set=c(3L, 2L, 2L))
)

test_that("download_OGRDB_germline_sequences()", {
    for (test_germline_set in .TEST_GERMLINE_SETS) {
        organism          <- test_germline_set$organism
        germline_sets     <- test_germline_set$germline_sets
        destfiles_per_set <- test_germline_set$destfiles_per_set

        ## --- download ungapped sequences ---

        dir.create(fasta_dir <- tempfile("fasta_dir_"))
        filenames <- download_OGRDB_germline_sequences(organism, germline_sets,
                                    gapped=FALSE, destdir=fasta_dir)
        expected_filenames <- paste0(
            head(igblastr:::IG_GROUPS, n=sum(destfiles_per_set)),
            ".fasta")
        expect_identical(unname(filenames), expected_filenames)
        expect_identical(names(filenames),
                         rep(names(germline_sets), destfiles_per_set))
        expect_true(setequal(list.files(fasta_dir), filenames))
        alleles <- readDNAStringSet(list.files(fasta_dir, full.names=TRUE))
        allele_names <- names(alleles)
        expect_true(is.character(allele_names))
        expect_equal(anyDuplicated(allele_names), 0)
        expect_false(anyNA(allele_names))

        ## --- check ungapped sequences against "coding sequences" ---
        ## ---    stored in the corresponding AIRR-C JSON files    ---

        dir.create(json_dir <- tempfile("json_dir_"))
        download_OGRDB_germline_json(organism, germline_sets,
                                     destdir=json_dir)
        json_paths <- list.files(json_dir, full.names=TRUE)
        df <- do.call(rbind, lapply(json_paths,
                                    extract_alleles_from_ogrdb_json))
        expect_equal(length(alleles), nrow(df))
        expect_true(setequal(allele_names, df$allele_name))
        m <- match(allele_names, df$allele_name)
        df <- S4Vectors:::extract_data_frame_rows(df, m)
        expected_alleles <- setNames(df$coding_sequence, df$allele_name)
        expect_identical(as.character(alleles), expected_alleles)

        ## --- download gapped sequences ---

        dir.create(fasta_dir2 <- tempfile("fasta_dir2_"))
        filenames2 <- download_OGRDB_germline_sequences(organism, germline_sets,
                                    gapped=TRUE, destdir=fasta_dir2)
        expect_identical(filenames2, filenames)
        expect_true(setequal(list.files(fasta_dir2), filenames2))
        alleles2 <- readDNAStringSet(list.files(fasta_dir2, full.names=TRUE))
        expect_identical(names(alleles2), allele_names)
        expect_true(all(width(alleles2) >= width(alleles)))
        expect_identical(width(alleles2) > width(alleles),
                         df[ , "sequence_type"] == "V")
        expect_identical(as.character(remove_gaps(alleles2)),
                         as.character(alleles))

        ## --- clean up ---

        unlink(fasta_dir, recursive=TRUE)
        unlink(json_dir, recursive=TRUE)
        unlink(fasta_dir2, recursive=TRUE)
    }
})


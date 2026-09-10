### =========================================================================
### install_OGRDB_germline_db()
### -------------------------------------------------------------------------


.collect_json_files_to_process <- function(json_dir, V_or_J)
{
    if (!isSingleNonWhiteString(json_dir))
        stop(wmsg("'json_dir' must be a single (non-empty) string"))
    if (!dir.exists(json_dir))
        stop(wmsg(json_dir, ": directory not found"))
    pattern <- sprintf("^(IG[HKL]|IG[HKL].*%s.*)\\.json$", V_or_J)
    json_files <- list.files(json_dir, pattern=pattern)
    if (length(json_files) == 0L)
        stop(wmsg(json_dir, ": directory contains no JSON files from ",
                  "which to extract auxiliary data"))
    json_files
}

.make_output_filenames_from_json_filenames <-
    function(json_filenames, output_suffix)
{
    stopifnot(is.character(json_filenames))
    if (!isSingleString(output_suffix))
        stop(wmsg("'output_suffix' must be a single string"))
    loci <- substr(json_filenames, 1L, 3L)
    stopifnot(all(loci %in% IG_LOCI))
    ## We don't allow more than one input JSON file per locus.
    if (anyDuplicated(loci)) {
        in1string <- paste0(json_filenames, collapse=", ")
        stop(wmsg("'json_dir' contains more than one ",
                  "input JSON file per locus: ", in1string))
    }
    add_suffix(loci, output_suffix)
}

.handle_existing_output_files <- function(output_files, overwrite, what)
{
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    existing_output_files <- output_files[file.exists(output_files)]
    if (length(existing_output_files) == 0L)
        return()
    if (overwrite) {
        unlink(existing_output_files)
        return()
    }
    s <- if (length(existing_output_files) >= 2L) "s" else ""
    in1string <- paste0(basename(existing_output_files), collapse=", ")
    msg1 <- c("'destdir' already contains the following ",
              "output ", what, " file", s, ": ", in1string)
    msg2 <- c("Use 'overwrite=TRUE' to overwrite them, or choose ",
              "another destination directory.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_intdata_files_from_ogrdb_jsons()
###

### 'json_dir' is assumed to contain JSON files downloaded with
### download_OGRDB_germline_json(). Only JSON files with names of the
### form <IG-locus>.json or <IG-locus>*V*.json are being considered for
### intdata extraction.
### Returns the names of the intdata files written to 'destdir' in a
### character vector. Note that this should **always** be a subset
### of 'c("IGH<output_suffix>", "IGK<output_suffix>", "IGL<output_suffix>")'.
make_intdata_files_from_ogrdb_jsons <-
    function(json_dir, destdir=".", output_suffix="V.ndm.imgt", overwrite=FALSE)
{
    json_files <- .collect_json_files_to_process(json_dir, "V")
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg(destdir, ": directory not found"))

    intdata_files <- .make_output_filenames_from_json_filenames(json_files,
                                                                output_suffix)

    loci <- substr(json_files, 1L, 3L)  # same as substr(intdata_files, 1L, 3L)

    intdata_files <- file.path(destdir, intdata_files)
    .handle_existing_output_files(intdata_files, overwrite, "intdata")

    json_files <- file.path(json_dir, json_files)
    ok <- vapply(seq_along(json_files),
        function(i) {
            intdata <- extract_intdata_from_ogrdb_json(json_files[[i]])
            if (nrow(intdata) == 0L)
                return(FALSE)
            write_ndm_data(intdata, intdata_files[[i]])
            TRUE
        }, logical(1))
    basename(intdata_files[ok])
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_auxdata_files_from_ogrdb_jsons()
###

.normarg_alleles_with_neg_cdr3_end <- function(alleles_with_neg_cdr3_end, loci)
{
    if (is.null(alleles_with_neg_cdr3_end))
        return(NULL)
    if (!is.list(alleles_with_neg_cdr3_end) ||
        length(alleles_with_neg_cdr3_end) == 0L)
        stop(wmsg("'alleles_with_neg_cdr3_end' must be NULL ",
                  "or a non-empty list"))
    nms <- names(alleles_with_neg_cdr3_end)
    if (is.null(nms) || anyDuplicated(nms))
        stop(wmsg("'alleles_with_neg_cdr3_end' must be supplied ",
                  "as a named list with unique names on it"))
    if (!all(nms %in% loci)) {
        in1string <- paste0('"', loci, '"', collapse=",")
        stop(wmsg("the names on 'alleles_with_neg_cdr3_end' must ",
                  "be a subset of 'c(", in1string, ")'"))
    }
    ok <- vapply(alleles_with_neg_cdr3_end, is.character, logical(1))
    if (!all(ok))
        stop(wmsg("the list elements in 'alleles_with_neg_cdr3_end' ",
                  "must be character vectors"))
    has_NAs <- vapply(alleles_with_neg_cdr3_end, anyNA, logical(1))
    if (any(has_NAs))
        stop(wmsg("the list elements in 'alleles_with_neg_cdr3_end' ",
                  "cannot contain NAs"))
    alleles_with_neg_cdr3_end
}

.stop_on_bad_alleles_not_as_expected <-
    function(found_bad_alleles, expected_bad_alleles)
{
    msg1 <- "Set of J alleles with negative cdr3 end is not as expected:"
    expected_in1string <- paste0(expected_bad_alleles, collapse=", ")
    found_in1string <- paste0(found_bad_alleles, collapse=", ")
    stop(wmsg2(msg1),
         "\n  - expected: ", wmsg2(expected_in1string, margin=14),
         "\n  -    found: ", wmsg2(found_in1string, margin=14))
}

.extract_and_fix_auxdata_from_ogrdb_json <- function(json_file, a_w_n_c_e=NULL)
{
    if (is.null(a_w_n_c_e))
        return(extract_auxdata_from_ogrdb_json(json_file))
    stopifnot(is.character(a_w_n_c_e), !anyNA(a_w_n_c_e))
    auxdata <- suppressWarnings(extract_auxdata_from_ogrdb_json(json_file))
    bad_idx <- which(auxdata[ , "cdr3_end"] < 0L)
    bad_alleles <- auxdata[bad_idx, "allele_name"]
    if (!setequal(bad_alleles, a_w_n_c_e))
        .stop_on_bad_alleles_not_as_expected(bad_alleles, a_w_n_c_e)
    auxdata[bad_idx, "cdr3_end"] <- NA_integer_
    auxdata
}

### 'json_dir' is assumed to contain JSON files downloaded with
### download_OGRDB_germline_json(). Only JSON files with names of the
### form <IG-locus>.json or <IG-locus>*J*.json are being considered for
### auxdata extraction.
### Returns the names of the auxdata files written to 'destdir' in a
### character vector. Note that this should **always** be a subset
### of 'c("IGH<output_suffix>", "IGK<output_suffix>", "IGL<output_suffix>")'.
make_auxdata_files_from_ogrdb_jsons <-
    function(json_dir, destdir=".", output_suffix="J_gl.aux",
             alleles_with_neg_cdr3_end=NULL, overwrite=FALSE)
{
    json_files <- .collect_json_files_to_process(json_dir, "J")
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg(destdir, ": directory not found"))

    auxdata_files <- .make_output_filenames_from_json_filenames(json_files,
                                                                output_suffix)

    loci <- substr(json_files, 1L, 3L)  # same as substr(auxdata_files, 1L, 3L)
    alleles_with_neg_cdr3_end <-
        .normarg_alleles_with_neg_cdr3_end(alleles_with_neg_cdr3_end, loci)

    auxdata_files <- file.path(destdir, auxdata_files)
    .handle_existing_output_files(auxdata_files, overwrite, "auxdata")

    json_files <- file.path(json_dir, json_files)
    ok <- vapply(seq_along(json_files),
        function(i) {
            a_w_n_c_e <- alleles_with_neg_cdr3_end[[loci[[i]]]]
            auxdata <- .extract_and_fix_auxdata_from_ogrdb_json(json_files[[i]],
                                                                a_w_n_c_e)
            if (nrow(auxdata) == 0L)
                return(FALSE)
            write_auxdata(auxdata, auxdata_files[[i]])
            TRUE
        }, logical(1))
    basename(auxdata_files[ok])
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_OGRDB_germline_db()
###

### Fetch all required FASTA and JSON files from OGRDB.
### Returns vector of corresponding loci.
.download_OGRDB_files_to_tmp_store <- function(organism, germline_sets,
                                               source_set, tmp_store, ...)
{
    fasta_files <- download_OGRDB_germline_sequences(organism, germline_sets,
                                      source_set=source_set,
                                      destdir=tmp_store, ...)
    loci <- unique(substr(fasta_files, 1L, 3L))
    json_files <- download_OGRDB_germline_json(organism, germline_sets,
                                      source_set=source_set,
                                      destdir=tmp_store, ...)
    stopifnot(setequal(unique(substr(json_files, 1L, 3L)), loci))
    IG_LOCI[IG_LOCI %in% loci]  # return loci in canonical order
}

.check_db_name_suffix <- function(suffix)
{
    if (!isSingleString(suffix))
        stop(wmsg("'suffix' must be a single (non-empty) string"))
    if (has_whitespace(suffix))
        stop(wmsg("'suffix' cannot contain whitespace characters"))
    if (grepl("[/\\]", suffix))
        stop(wmsg("'suffix' cannot contain slahes, forward (/) or ",
                  "backward (\\)"))
}

### We assume that the longest string contains the most specific strain
### e.g. "C57BL/6J". All the other strings must be prefixes of the longest
### string e.g. "" or "C57BL/6". If that's not the case, then it means that
### the user picked up incompatible mouse germline sets e.g. they mixed
### C57BL/6 sets with CAST/EiJ sets.
.most_specific_strain <- function(strains)
{
    stopifnot(is.character(strains), length(strains) != 0L)
    nc <- nchar(strains)
    ans <- unique(strains[nc == max(nc)])
    if (length(ans) >= 2L || !all(startsWith(ans, strains)))
        stop(wmsg("cannot mix germline sets from different mouse strains"))
    ans
}

.infer_mouse_strain_from_OGRDB_set_names <- function(set_names)
{
    strains <- extract_mouse_strains_from_OGRDB_set_names(set_names)
    strain <- .most_specific_strain(strains)
    chartr("/", "_", strain)
}

.form_OGRDB_germline_db_name <- function(organism, set_names,
                                         loci, suffix, source_set)
{
    .check_db_name_suffix(suffix)
    if (organism == "Mus musculus") {
        strain <- .infer_mouse_strain_from_OGRDB_set_names(set_names)
        if (nzchar(strain))
            organism <- paste0(organism, ".", strain)
    }
    organism <- chartr(" ", "_", organism)
    db_name <- sprintf("OGRDB.%s.%s", organism, paste(loci, collapse="+"))
    if (nzchar(suffix))
        db_name <- paste0(db_name, ".", suffix)
    if (source_set)
        db_name <- paste0(db_name, ".src")
    db_name
}

.extract_intdata_from_ogrdb_jsons <- function(tmp_store)
{
    make_intdata_files_from_ogrdb_jsons(tmp_store, destdir=tmp_store,
                                        output_suffix="V.ndm.imgt")
    intdata_files <- list.files(tmp_store, pattern="V.ndm.imgt$",
                                full.names=TRUE)
    do.call(rbind, lapply(intdata_files, read_ndm_data))
}

.extract_auxdata_from_ogrdb_jsons <- function(tmp_store)
{
    make_auxdata_files_from_ogrdb_jsons(tmp_store, destdir=tmp_store,
                                        output_suffix="J_gl.aux")
    auxdata_files <- list.files(tmp_store, pattern="J_gl.aux$",
                                full.names=TRUE)
    do.call(rbind, lapply(auxdata_files, read_auxdata))
}

install_OGRDB_germline_db <- function(organism, germline_sets,
                                      suffix="", source_set=FALSE,
                                      overwrite=FALSE, verbose=FALSE, ...)
{
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    tmp_store <- tempfile()
    dir.create(tmp_store)
    on.exit(nuke_file(tmp_store))

    loci <- .download_OGRDB_files_to_tmp_store(organism, germline_sets,
                                               source_set, tmp_store, ...)

    ## Prepare 'db_name'.
    db_name <- .form_OGRDB_germline_db_name(organism, names(germline_sets),
                                            loci, suffix, source_set)

    ## Prepare 'intdata' and 'auxdata'.
    intdata <- .extract_intdata_from_ogrdb_jsons(tmp_store)
    auxdata <- .extract_auxdata_from_ogrdb_jsons(tmp_store)

    ## Create and install germline db.
    install_dir <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    if.exists <- if (overwrite) "overwrite" else "error"
    install_germline_db(install_dir, db_name, tmp_store, loci,
                        gapped=TRUE, intdata=intdata, auxdata=auxdata,
                        if.exists=if.exists, verbose=verbose,
                        cheer.if.success=TRUE)

    invisible(db_name)
}


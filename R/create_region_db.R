### =========================================================================
### create_region_db() and family
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.checkarg_fasta_files <- function(fasta_files)
{
    if (!is.character(fasta_files) || length(fasta_files) == 0L)
        stop(wmsg("'fasta_files' must be a non-empty character vector"))
    if (anyNA(fasta_files) || anyDuplicated(fasta_files))
        stop(wmsg("'fasta_files' cannot contain NAs or duplicates"))
    if (any(is_white_str(fasta_files)))
        stop(wmsg("'fasta_files' cannot contain empty or white strings"))
}

.checkarg_destdir <- function(destdir)
{
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
}

.get_db_final_fasta_path <- function(destdir, region_type)
{
    .checkarg_destdir(destdir)
    file.path(destdir, paste0(region_type, ".fasta"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .load_allele_set()
###

.infer_loci_from_filenames <- function(fasta_files, region_type)
{
    .checkarg_fasta_files(fasta_files)
    stopifnot(isSingleNonWhiteString(region_type), nchar(region_type) == 1L)
    ## If 'fasta_files' has names on it, we use them as the actual file names.
    filenames <- names(fasta_files)
    if (is.null(filenames))
        filenames <- basename(fasta_files)
    ## The strings in 'filenames' must be of the form <locus>V.fasta
    ## or <locus>J.fasta where locus is a 3-letter string. Note that we don't
    ## really care if <locus> is a valid locus or not (see IG_LOCI and TR_LOCI
    ## in R/loci-utils.R for the 7 valid loci), as long as it's a 3-letter
    ## string, so we don't bother to check.
    suffixes <- substr(filenames, 4L, nchar(filenames))
    expected_suffix <- paste0(region_type, ".fasta")
    if (!all(suffixes == expected_suffix)) {
        in1string <- paste(filenames, collapse=", ")
        stop(wmsg("The names of the FASTA files must be of the form ",
                  "<locus>", expected_suffix, ", where <locus> is a ",
                  "3-letter string. Got: ", in1string))
    }
    substr(filenames, 1L, 3L)
}

.load_allele_set <- function(fasta_files, region_type,
                             excluded_alleles=character(0), verbose=FALSE)
{
    loci <- .infer_loci_from_filenames(fasta_files, region_type)
    stopifnot(is.character(excluded_alleles), isTRUEorFALSE(verbose))

    if (verbose) {
        in1string <- paste(basename(fasta_files), collapse=", ")
        msg <- c("Input file(s): ", in1string)
        message("  o ", wmsg(msg, margin=4L), " ", appendLF=FALSE)
    }

    allele_sets <- lapply(seq_along(fasta_files),
        function(i) {
            allele_set <- readDNAStringSet(fasta_files[[i]])
            mcols(allele_set)$locus <- loci[[i]]
            allele_set
        }
    )
    allele_set <- do.call(c, allele_sets)
    if (length(allele_set) == 0L) {
        in1string <- paste(basename(fasta_files), collapse=", ")
        stop(wmsg("no alleles found in FASTA files: ", in1string))
    }

    if (verbose)
        message("(", length(allele_set), " alleles)\n")

    if (length(excluded_alleles) != 0L) {
        allele_names <- clean_imgt_fasta_headers(names(allele_set))
        excluded_idx <- which(allele_names %in% excluded_alleles)
        excluded_idx_len <- length(excluded_idx)
        if (excluded_idx_len != 0L) {
            if (verbose) {
                in1string <- paste(allele_names[excluded_idx], collapse=", ")
                msg <- c("Removed ", excluded_idx_len, " explicitly ",
                         "excluded allele(s) from initial set of ",
                         length(allele_set), " alleles: ", in1string, ".")
                message("  o ", wmsg(msg, margin=4L), "\n")
            }
            allele_set <- allele_set[-excluded_idx]
        }
    }
    allele_set
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .init_region_db()
###

.region_db_already_exists <- function(destdir, region_type)
{
    original_fasta_dir <- get_db_original_fasta_dir(destdir, region_type)
    final_fasta <- .get_db_final_fasta_path(destdir, region_type)
    file.exists(original_fasta_dir) || file.exists(final_fasta)
}

.stop_on_existing_region_db <- function(destdir, region_type)
{
    msg1 <- c("There seems to be already a ", region_type, "-region ",
              "database in ", destdir)
    msg2 <- c("Use 'overwrite=TRUE' to overwrite or choose another ",
              "destination directory.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2))
}

.nuke_existing_region_db <- function(destdir, region_type)
{
    pattern <- paste0("^", region_type, "\\.fasta$")
    clean_blastdbs(destdir, pattern)
    original_fasta_dir <- get_db_original_fasta_dir(destdir, region_type)
    final_fasta <- .get_db_final_fasta_path(destdir, region_type)
    nuke_file(original_fasta_dir)
    nuke_file(final_fasta)
}

.init_region_db <- function(fasta_files, destdir, region_type,
                            excluded_alleles=character(0),
                            overwrite=FALSE, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
    if (verbose)
        message("Creating ", region_type, "-region db ...\n")

    .checkarg_fasta_files(fasta_files)
    .checkarg_destdir(destdir)
    if (!is.character(excluded_alleles))
        stop(wmsg("'excluded_alleles' must be a character vector"))
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (.region_db_already_exists(destdir, region_type)) {
        if (!overwrite)
            .stop_on_existing_region_db(destdir, region_type)
        .nuke_existing_region_db(destdir, region_type)
    }

    .load_allele_set(fasta_files, region_type,
                     excluded_alleles=excluded_alleles, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .write_db_fasta_files()
###

### Write the clean allele set to the final FASTA file.
.write_db_final_fasta <- function(cleaned_allele_set, destdir, region_type,
                                  verbose=FALSE)
{
    final_fasta <- .get_db_final_fasta_path(destdir, region_type)
    if (verbose) {
        msg <- c("Writing final FASTA file (", basename(final_fasta), ") ",
                 "to the db")
        message("  o ", wmsg(msg, margin=4L), " ... ", appendLF=FALSE)
    }
    writeXStringSet(cleaned_allele_set, final_fasta)
    if (verbose)
        message("ok.\n")
}

### Returns a character vector containing the allele names in the
### input FASTA file ('infasta').
.copy_fasta_file <- function(infasta, outfasta, DORsummary)
{
    stopifnot(isSingleNonWhiteString(infasta),
              !dir.exists(infasta), file.exists(infasta),
              isSingleNonWhiteString(outfasta),
              !dir.exists(outfasta), !file.exists(outfasta),
              is.data.frame(DORsummary))
    headers <- names(fasta.seqlengths(infasta))
    allele_names <- clean_imgt_fasta_headers(headers)

    ## We use the "Drop Or Rename summary" ('DORsummary') to identify the
    ## FASTA records to drop.
    ## There are two reasons for dropping a FASTA record:
    ##   1. The record is for an allele that belongs to the exclusion list
    ##      i.e. the allele got excluded early by .init_region_db() based
    ##      on the list of allele names passed to its 'excluded_alleles'
    ##      argument (see function .init_region_db() above in this file).
    ##      This means that the allele didn't make it to the 'allele_set'
    ##      object (DNAStringSet) returned by .init_region_db(), so is not
    ##      in 'DORsummary$allele_name' either.
    ##   2. The record is for an allele that made it to the 'allele_set' object
    ##      returned by .init_region_db(), but the allele was later marked
    ##      as "repeated" by clean_allele_set() or clean_[VJ]_allele_set().
    ##      In this case the allele name is in 'DORsummary$allele_name' but
    ##      with a corresponding suffix in 'DORsummary$suffix' set to NA.
    drop_me <- !(allele_names %in% DORsummary[ , "allele_name"])
    stopifnot(identical(allele_names[!drop_me], DORsummary[ , "allele_name"]))
    drop_me[!drop_me] <- is.na(DORsummary[ , "suffix"])

    drop_idx <- which(drop_me)
    if (length(drop_idx) == 0L) {
        stopifnot(file.copy(infasta, outfasta))
    } else {
        ## Drop FASTA records. The only reason we use a readLines/writeLines
        ## approach instead of a readBStringSet/writeBStringSet approach is
        ## because we want to copy the FASTA records in their original form
        ## (the readBStringSet/writeBStringSet approach would reformat the
        ## multi-line sequences into fixed-width lines). This would allow us
        ## to use Unix command diff to compare the input and output files, if
        ## we wanted to, and get a clean/minimalist output that shows what was
        ## dropped. The drawback of the readLines/writeLines approach is that
        ## it's quite inefficient on big files (although the FASTA files we
        ## use to create a region db are typically small).
        lines <- readLines(infasta)
        record_starts <- grep("^>", lines)
        stopifnot(identical(lines[record_starts], paste0(">", headers)))
        record_ends <- c(tail(record_starts - 1L, n=-1L), length(lines))
        record_bounds <- IRanges(record_starts, record_ends)
        lines <- extractROWS(lines, record_bounds[-drop_idx])
        writeLines(lines, outfasta)
    }
    allele_names
}

### Returns a character vector containing the allele names in the
### input FASTA files ('fasta_files').
.copy_fasta_files_to_db_original_fasta_dir <-
    function(fasta_files, original_fasta_dir, DORsummary)
{
    .checkarg_fasta_files(fasta_files)
    out_filenames <- names(fasta_files)
    if (is.null(out_filenames))
        out_filenames <- basename(fasta_files)
    stopifnot(all(has_suffix(out_filenames, ".fasta")),
              isSingleNonWhiteString(original_fasta_dir),
              dir.exists(original_fasta_dir),
              is.data.frame(DORsummary))

    ## Split 'DORsummary'.
    locus <- DORsummary[ , "locus"]
    check_locus(locus, "DORsummary$locus")
    f <- factor(locus, levels=substr(out_filenames, 1L, 3L))
    stopifnot(!anyNA(f))
    DORsummaries <- split(DORsummary, f)

    all_allele_names <- lapply(seq_along(fasta_files),
        function(i) {
            infasta <- fasta_files[[i]]
            outfasta <- file.path(original_fasta_dir, out_filenames[[i]])
            .copy_fasta_file(infasta, outfasta, DORsummaries[[i]])
        }
    )
    unlist(all_allele_names, use.names=FALSE)
}

### Create "original fasta" subdir and copy fasta files to it.
### 'excluded_alleles' used for sanity check only.
.make_db_original_fasta_dir <-
    function(fasta_files, destdir, region_type,
             DORsummary, excluded_alleles=character(0))
{
    original_fasta_dir <- get_db_original_fasta_dir(destdir, region_type)
    stopifnot(!file.exists(original_fasta_dir), dir.create(original_fasta_dir))
    allele_names <- .copy_fasta_files_to_db_original_fasta_dir(fasta_files,
                                            original_fasta_dir, DORsummary)
    ## Sanity checks.
    DORsummary_allele_name <- DORsummary[ , "allele_name"]
    stopifnot(all(DORsummary_allele_name %in% allele_names))
    was_excluded <- !(allele_names %in% DORsummary_allele_name)
    stopifnot(identical(allele_names[!was_excluded], DORsummary_allele_name),
              all(allele_names[was_excluded] %in% excluded_alleles))
}

.write_db_fasta_files <-
    function(cleaned_allele_set, destdir, region_type,
             fasta_files, excluded_alleles=character(0),
             verbose=FALSE)
{
    DORsummary <- metadata(cleaned_allele_set)$DORsummary
    stopifnot(is.data.frame(DORsummary))
    suffix <- DORsummary[ , "suffix"]
    expected_allele_names <- add_suffix(DORsummary[ , "allele_name"], suffix)
    expected_allele_names <- expected_allele_names[!is.na(suffix)]
    stopifnot(identical(names(cleaned_allele_set), expected_allele_names))

    .write_db_final_fasta(cleaned_allele_set, destdir, region_type,
                          verbose=verbose)
    .make_db_original_fasta_dir(fasta_files, destdir, region_type,
                                DORsummary, excluded_alleles=excluded_alleles)
    if (verbose)
        message("... done. (Final number of ", region_type, " alleles ",
                "in db = ", length(cleaned_allele_set), ")\n")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_region_db()
###
### Generic creation of a "region db" (V-, D-, J-, or C-region) from a
### collection of FASTA files (typically obtained from IMGT or
### AIRR-community/OGRDB) for a given organism.
###
### See procedure described at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### for how to create a germline or C-region db from the FASTA files
### available at IMGT.
### This is a 3-step procedure: (1) combine, (2) clean, (3) compile.
### create_region_db() only performs steps (1) and (2). Compilation (with
### makeblastdb) will happen at a latter time.
###
### See create_V_region_db() and create_J_region_db() below for more
### specialized versions.

### 'destdir' must be the path to a writable directory that already exists!
### See clean_allele_set() in R/clean_allele_set.R for the role of
### the 'disambiguate.allele.names' argument.
### The following subdirectory and file will be added to 'destdir':
### - D_original_fasta/: subdirectory containing the input FASTA files
###       corresponding to the D regions, one FASTA file per region;
### - D.fasta: final FASTA file containing the clean allele set.
create_region_db <- function(destdir, fasta_files,
                             region_type=VDJC_REGION_TYPES,
                             excluded_alleles=character(0),
                             disambiguate.allele.names=FALSE,
                             overwrite=FALSE, verbose=FALSE)
{
    region_type <- match.arg(region_type)
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))

    ## (1) Load and combine.
    allele_set <- .init_region_db(fasta_files, destdir, region_type,
                                  excluded_alleles=excluded_alleles,
                                  overwrite=overwrite, verbose=verbose)

    ## (2) Clean.
    cleaned_allele_set <-
        clean_allele_set(allele_set,
                         disambiguate.allele.names=disambiguate.allele.names,
                         verbose=verbose)

    ## Write the db.
    .write_db_fasta_files(cleaned_allele_set, destdir, region_type,
                          fasta_files, excluded_alleles=excluded_alleles,
                          verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_V_region_db()
###

.add_auto_intdata_to_db <- function(cleaned_allele_set, destdir, verbose=FALSE)
{
    if (verbose) {
        msg <- "Adding the igblastr-generated intdata to the db"
        message("  o ", wmsg(msg, margin=4L), " ... ", appendLF=FALSE)
    }
    stopifnot(is(cleaned_allele_set, "DNAStringSet"))
    allele_names <- names(cleaned_allele_set)
    stopifnot(!is.null(allele_names))
    cleaned_allele_set_mcols <- mcols(cleaned_allele_set, use.names=FALSE)
    stopifnot(
        is(cleaned_allele_set_mcols, "DataFrame"),
        identical(allele_names, cleaned_allele_set_mcols[ , "allele_name"])
    )
    ndm_data <- cleaned_allele_set_mcols[ , names(NDM_DATA_COL2CLASS)]
    write_intdata_to_db(as.data.frame(ndm_data), destdir)
    if (verbose)
        message("ok.\n")
}

### A specialized version of create_region_db() for V alleles.
### Set 'gapped' to TRUE if the supplied sequences are gapped V allele
### sequences (note that only V allele sequences are allowed to have gaps).
### Set 'auto.intdata' to TRUE if the supplied sequences are gapped V allele
### sequences **and** the associated internal data should be computed and
### added to the db.
create_V_region_db <- function(destdir, fasta_files,
                               gapped=FALSE, auto.intdata=FALSE,
                               fwrcdr_widths=IMGT_FWRCDR_WIDTHS,
                               disambiguate.allele.names=FALSE,
                               overwrite=FALSE, verbose=FALSE)
{
    if (!isTRUEorFALSE(gapped))
        stop(wmsg("'gapped' must be TRUE or FALSE"))
    checkarg_auto.intdata(auto.intdata, gapped)
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))

    ## (1) Load and combine.
    ## Note that 'allele_set' is loaded with the locus info as a metadata
    ## column on it. clean_V_allele_set() will require and use this metadata
    ## column when 'auto.intdata' is TRUE.
    allele_set <- .init_region_db(fasta_files, destdir, "V",
                                  overwrite=overwrite, verbose=verbose)

    ## (2) Clean.
    cleaned_allele_set <-
        clean_V_allele_set(allele_set,
                           gapped=gapped, auto.intdata=auto.intdata,
                           fwrcdr_widths=fwrcdr_widths,
                           disambiguate.allele.names=disambiguate.allele.names,
                           verbose=verbose)

    ## Write the db.
    if (auto.intdata)
        .add_auto_intdata_to_db(cleaned_allele_set, destdir, verbose=verbose)
    .write_db_fasta_files(cleaned_allele_set, destdir, "V",
                          fasta_files, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_J_region_db()
###

### Adds the igblastr-generated auxdata to the db ONLY if
### its "coding_frame_start" column contains no NAs.
.add_auto_auxdata_to_db <- function(auxdata, destdir, verbose=FALSE)
{
    bad_alleles <- J_alleles_with_missing_coding_frame_start(auxdata)
    if (length(bad_alleles) != 0L) {
        if (verbose) {
            in1string <- paste(bad_alleles, collapse=", ")
            message(wmsg2("--> compute_auxdata() could not determine ",
                          "the \"coding frame start\" for J allele(s): ",
                          in1string, margin=4L))
            message(wmsg("--> NOT adding the igblastr-generated auxdata ",
                         "to the db."), "\n")
        }
        return()
    }

    if (verbose) {
        msg <- "Adding the igblastr-generated auxdata to the db"
        message("  o ", wmsg(msg, margin=4L), " ... ", appendLF=FALSE)
    }

    auxdata <- auxdata[ , names(AUXDATA_COL2CLASS)]

    ## write_auxdata_to_db() will reject a data.frame with negative
    ## values in the "cdr3_end" column, so we replace them with NAs.
    bad_idx <- which(auxdata[ , "cdr3_end"] < 0L)
    auxdata[bad_idx, "cdr3_end"] <- NA_integer_
    write_auxdata_to_db(auxdata, destdir)

    if (verbose)
        message("ok.\n")
}

### A specialized version of create_region_db() for J alleles.
### Set 'auto.auxdata' to TRUE if the auxiliary data associated with the
### J alleles should be computed and added to the db.
create_J_region_db <- function(destdir, fasta_files, imgt.fasta.headers=FALSE,
                               excluded_alleles=character(0),
                               auto.auxdata=FALSE, ref_auxdata=NULL,
                               disambiguate.allele.names=FALSE,
                               overwrite=FALSE, verbose=FALSE)
{
    if (!isTRUEorFALSE(imgt.fasta.headers))
        stop(wmsg("'imgt.fasta.headers' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(auto.auxdata))
        stop(wmsg("'auto.auxdata' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))

    ## (1) Load and combine.
    ## Note that 'allele_set' is loaded with the locus info as a metadata
    ## column on it. clean_J_allele_set() will require and use this metadata
    ## column when 'auto.auxdata' is TRUE.
    allele_set <- .init_region_db(fasta_files, destdir, "J",
                                  excluded_alleles=excluded_alleles,
                                  overwrite=overwrite, verbose=verbose)

    ## (2) Clean.
    cleaned_allele_set <-
        clean_J_allele_set(allele_set, imgt.fasta.headers=imgt.fasta.headers,
                           auto.auxdata=auto.auxdata, ref_auxdata=ref_auxdata,
                           disambiguate.allele.names=disambiguate.allele.names,
                           verbose=verbose)

    ## Write the db.
    if (auto.auxdata) {
        allele_names <- names(cleaned_allele_set)
        auxdata <- mcols(cleaned_allele_set, use.names=FALSE)
        stopifnot(identical(allele_names, auxdata[ , "allele_name"]))
        auxdata <- as.data.frame(auxdata)
        .add_auto_auxdata_to_db(auxdata, destdir, verbose=verbose)
    }
    .write_db_fasta_files(cleaned_allele_set, destdir, "J",
                          fasta_files, excluded_alleles=excluded_alleles,
                          verbose=verbose)
}


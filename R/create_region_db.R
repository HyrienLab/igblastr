### =========================================================================
### create_region_db() and related
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
}

.checkarg_destdir <- function(destdir)
{
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
}

.get_final_fasta_path <- function(destdir, region_type)
{
    .checkarg_destdir(destdir)
    file.path(destdir, paste0(region_type, ".fasta"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers for create_region_db() and related
###

.region_db_already_exists <- function(destdir, region_type)
{
    original_fasta_dir <- get_db_original_fasta_dir(destdir, region_type)
    final_fasta <- .get_final_fasta_path(destdir, region_type)
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
    final_fasta <- .get_final_fasta_path(destdir, region_type)
    nuke_file(original_fasta_dir)
    nuke_file(final_fasta)
}

### Creates "original fasta" subdir and copy fasta files to it.
### Returns the paths to the copied files.
.init_region_db <- function(fasta_files, destdir, region_type,
                            overwrite=FALSE, verbose=FALSE)
{
    if (verbose)
        message("Creating ", region_type, "-region db ...\n")
    .checkarg_fasta_files(fasta_files)
    .checkarg_destdir(destdir)
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
    if (.region_db_already_exists(destdir, region_type)) {
        if (!overwrite)
            .stop_on_existing_region_db(destdir, region_type)
        .nuke_existing_region_db(destdir, region_type)
    }
    original_fasta_dir <- get_db_original_fasta_dir(destdir, region_type)
    stopifnot(dir.create(original_fasta_dir))
    copy_files_to_dir(fasta_files, original_fasta_dir)
    original_fasta_files <- list_fasta_files(original_fasta_dir)
    if (verbose) {
        in1string <- paste(basename(original_fasta_files), collapse=", ")
        msg <- c("Input file(s): ", in1string)
        message("  o ", wmsg(msg, margin=4L), "\n")
    }
    original_fasta_files
}

### Write the clean allele set to the final FASTA file.
.dump_clean_allele_set <- function(allele_set, destdir, region_type,
                                   verbose=FALSE)
{
    final_fasta <- .get_final_fasta_path(destdir, region_type)
    writeXStringSet(allele_set, final_fasta)
    if (verbose)
        message("... done. (Final number of ", region_type, " alleles ",
                "in db = ", length(allele_set), ")\n")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .load_original_alleles_with_loci()
###

### The names of the supplied files must be of the form <locus>V.fasta
### or <locus>J.fasta where locus is a 3-letter string. Note that we don't
### really care if <locus> is a valid locus or not (see IG_LOCI and TR_LOCI
### in R/loci-utils.R for the 7 valid loci), as long as it's a 3-letter
### string, so we don't bother to check.
.infer_loci_from_filenames <- function(fasta_files, region_type)
{
    .checkarg_fasta_files(fasta_files)
    filenames <- basename(fasta_files)
    suffix <- substr(filenames, 4L, nchar(filenames))
    expected_suffix <- paste0(region_type, ".fasta")
    if (!all(suffix == expected_suffix)) {
        in1string <- paste(filenames, collapse=", ")
        stop(wmsg("The names of the FASTA files must be of the form ",
                  "<locus>", expected_suffix, ", where <locus> is a ",
                  "3-letter string. Got: ", in1string))
    }
    substr(filenames, 1L, 3L)
}

.load_original_alleles_with_loci <- function(fasta_files, region_type)
{
    loci <- .infer_loci_from_filenames(fasta_files, region_type)
    dna <- lapply(seq_along(fasta_files),
        function(i) {
            dna <- readDNAStringSet(fasta_files[[i]])
            mcols(dna)$locus <- loci[[i]]
            dna
        }
    )
    do.call(c, dna)
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
create_region_db <- function(fasta_files, destdir,
                             region_type=VDJC_REGION_TYPES,
                             disambiguate.allele.names=FALSE,
                             overwrite=FALSE, verbose=FALSE)
{
    region_type <- match.arg(region_type)
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))

    ## Create "original fasta" subdir and copy fasta files to it.
    original_fasta_files <- .init_region_db(fasta_files, destdir, region_type,
                                            overwrite=overwrite,
                                            verbose=verbose)

    ## (1) Load and combine.
    allele_set <- readDNAStringSet(original_fasta_files)

    ## (2) Clean.
    allele_set <- clean_allele_set(allele_set,
                        disambiguate.allele.names=disambiguate.allele.names,
                        verbose=verbose)

    ## Write the clean allele set to the final FASTA file.
    .dump_clean_allele_set(allele_set, destdir, region_type, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_V_region_db()
###

.add_intdata_to_db <- function(allele_set, destdir)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))
    allele_set_mcols <- mcols(allele_set, use.names=FALSE)
    stopifnot(is(allele_set_mcols, "DataFrame"),
              identical(allele_names, allele_set_mcols[ , "allele_name"]))

    intdata_path <- make_germline_db_intdata_path(destdir, for.aa=FALSE,
                                                  domain_system="imgt")
    intdata_dir <- dirname(intdata_path)
    stopifnot(!dir.exists(intdata_dir))

    ndm_data <- as.data.frame(allele_set_mcols[ , names(NDM_DATA_COL2CLASS)])
    check_ndm_data_col2class(ndm_data)
    dir.create(intdata_dir)
    write_ndm_data(ndm_data, intdata_path)
}

### A specialized version of create_region_db() for V alleles.
### Set 'gapped' to TRUE if the supplied sequences are gapped V allele
### sequences (note that only V allele sequences are allowed to have gaps).
### Set 'with.intdata' to TRUE if the supplied sequences are gapped V allele
### sequences **and** the associated internal data should be computed and
### added to the db.
create_V_region_db <- function(fasta_files, destdir,
                               gapped=FALSE, with.intdata=FALSE,
                               disambiguate.allele.names=FALSE,
                               overwrite=FALSE, verbose=FALSE)
{
    if (!isTRUEorFALSE(gapped))
        stop(wmsg("'gapped' must be TRUE or FALSE"))
    checkarg_with.intdata(with.intdata, gapped)
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))

    ## Create "original fasta" subdir and copy fasta files to it.
    original_fasta_files <- .init_region_db(fasta_files, destdir, "V",
                                            overwrite=overwrite,
                                            verbose=verbose)

    ## (1) Load and combine.
    if (with.intdata) {
        ## We use .load_original_alleles_with_loci() to get the locus info
        ## (returned as a metadata column on 'allele_set').
        ## clean_V_allele_set() requires and uses this metadata column
        ## when 'with.intdata' is TRUE.
        allele_set <-
            .load_original_alleles_with_loci(original_fasta_files, "V")
    } else {
        allele_set <- readDNAStringSet(original_fasta_files)
    }

    ## (2) Clean.
    allele_set <- clean_V_allele_set(allele_set,
                        gapped=gapped, with.intdata=with.intdata,
                        disambiguate.allele.names=disambiguate.allele.names,
                        verbose=verbose)

    ## Write the clean allele set to the final FASTA file.
    .dump_clean_allele_set(allele_set, destdir, "V", verbose=verbose)

    if (with.intdata) {
        if (verbose)
            message(wmsg("Adding the intdata to the db"), " ... ",
                    appendLF=FALSE)
        .add_intdata_to_db(allele_set, destdir)
        if (verbose)
            message("ok.\n")
    }

    allele_set
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_J_region_db()
###

### Not exported!
### Adds the auxdata to the db ONLY if its "coding_frame_start" column
### contains no NAs.
### Returns TRUE if that's actually the case (and therefore the auxdata
### was added), and FALSE otherwise.
add_computed_auxdata_to_db <- function(auxdata, destdir, verbose=FALSE)
{
    bad_alleles <- J_alleles_with_missing_coding_frame_start(auxdata)
    if (length(bad_alleles) != 0L) {
        if (verbose) {
            in1string <- paste(bad_alleles, collapse=", ")
            message(wmsg2("--> compute_auxdata() could not determine ",
                          "the \"coding frame start\" for J allele(s): ",
                          in1string, margin=4L))
            message(wmsg("--> NOT adding the auxdata to the db."), "\n")
        }
        return(FALSE)
    }

    if (verbose)
        message(wmsg("Adding the computed auxdata to the db"), " ... ",
                appendLF=FALSE)

    ## Even though write_auxdata() will call check_auxdata_col2class()
    ## internally, we prefer to fail **before** creating the 'auxdata_dir'
    ## folder.
    auxdata <- auxdata[ , names(AUXDATA_COL2CLASS)]
    check_auxdata_col2class(auxdata)

    auxdata_path <- make_germline_db_auxdata_path(destdir)
    auxdata_dir <- dirname(auxdata_path)
    stopifnot(!dir.exists(auxdata_dir))
    dir.create(auxdata_dir)

    write_auxdata(auxdata, auxdata_path)

    if (verbose)
        message("ok.\n")
    TRUE
}

### A specialized version of create_region_db() for J alleles.
### Set 'with.auxdata' to TRUE if the auxiliary data associated with the
### J alleles should be computed and added to the db.
create_J_region_db <- function(fasta_files, destdir,
                               excluded_J_alleles=NULL,
                               with.auxdata=FALSE, imgt.fasta=FALSE,
                               disambiguate.allele.names=FALSE,
                               overwrite=FALSE, verbose=FALSE)
{
    if (!isTRUEorFALSE(with.auxdata))
        stop(wmsg("'with.auxdata' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(imgt.fasta))
        stop(wmsg("'imgt.fasta' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))

    ## Create "original fasta" subdir and copy fasta files to it.
    original_fasta_files <- .init_region_db(fasta_files, destdir, "J",
                                            overwrite=overwrite,
                                            verbose=verbose)

    ## (1) Load and combine.
    if (with.auxdata) {
        ## We use .load_original_alleles_with_loci() to get the locus info
        ## (returned as a metadata column on 'allele_set').
        ## clean_J_allele_set() requires and uses this metadata column
        ## when 'with.auxdata' is TRUE.
        allele_set <-
            .load_original_alleles_with_loci(original_fasta_files, "J")
    } else {
        allele_set <- readDNAStringSet(original_fasta_files)
    }

    ## (2) Clean.
    allele_set <- clean_J_allele_set(allele_set,
                        excluded_J_alleles=excluded_J_alleles,
                        with.auxdata=with.auxdata, imgt.fasta=imgt.fasta,
                        disambiguate.allele.names=disambiguate.allele.names,
                        verbose=verbose)

    ## Write the clean allele set to the final FASTA file.
    .dump_clean_allele_set(allele_set, destdir, "J", verbose=verbose)

    if (with.auxdata) {
        allele_names <- names(allele_set)
        auxdata <- mcols(allele_set, use.names=FALSE)
        stopifnot(identical(allele_names, auxdata[ , "allele_name"]))
        auxdata <- as.data.frame(auxdata)
        add_computed_auxdata_to_db(auxdata, destdir, verbose=verbose)
    }

    allele_set
}


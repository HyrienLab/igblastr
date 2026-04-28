### =========================================================================
### Add novel gene alleles to a germline db
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .novel_alleles_as_temp_edited_fasta()
###

### 'novel_alleles' must be the path to a FASTA file (possibly
### gz-compressed), or a named DNAStringSet object.
### Returns **absolute** path to the **uncompressed** FASTA file
### containing the germline sequences.
.normarg_novel_alleles <- function(novel_alleles)
{
    if (isSingleNonWhiteString(novel_alleles))
        return(fasta_files_as_one_uncompressed_file(novel_alleles,
                                                    "novel_alleles"))
    if (is(novel_alleles, "DNAStringSet")) {
        if (length(novel_alleles) != 0L) {
            if (is.null(names(novel_alleles)))
                stop(wmsg("DNAStringSet object 'novel_alleles' ",
                          "must have names"))
            seqlens <- setNames(width(novel_alleles), names(novel_alleles))
            check_seqlens(seqlens, "query")
        }
        path <- tempfile("novel_alleles_", fileext=".fasta")
        writeXStringSet(novel_alleles, path)
        attr(path, "safe_to_remove") <- TRUE
        return(path)
    }
    stop(wmsg("'novel_alleles' must be a single (non-empty) string ",
              "that is the path to a FASTA file, or a named DNAStringSet ",
              "object"))
}

### Checks and returns the novel alleles in a temp FASTA file.
.novel_alleles_as_temp_edited_fasta <- function(novel_alleles)
{
    novel_fasta <- .normarg_novel_alleles(novel_alleles)
    if (isTRUE(attr(novel_fasta, "safe_to_remove")))
        on.exit(unlink(novel_fasta))
    edited_novel_fasta <- tempfile("edited_novel_fasta_", fileext = ".fasta")
    redit_imgt_file(novel_fasta, edited_novel_fasta)
    edited_novel_fasta
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .check_novel_allele_names()
###

### 'alleles' must be the path to a FASTA file that is assumed to have been
### processed with redit_imgt_file().
.extract_allele_names <- function(alleles, in_what)
{
    stopifnot(isSingleNonWhiteString(alleles))
    allele_names <- names(fasta.seqlengths(alleles))
    allele_names <- trimws2(sub("\\|.*$", "", allele_names))
    if (!all(nzchar(allele_names)))
        stop(wmsg("all the sequences in ", in_what, " must have a name"))
    dupidx <- which(duplicated(allele_names))
    if (length(dupidx) != 0L) {
        in1string <- paste(allele_names[dupidx], collapse=", ")
        stop(wmsg("the following allele names in ", in_what, " are ",
                  "duplicated: ", in1string))
    }
    allele_names
}

.check_novel_allele_names <- function(db_fasta_file, db_name, novel_alleles)
{
    stopifnot(isSingleNonWhiteString(db_fasta_file),
              isSingleNonWhiteString(db_name),
              isSingleNonWhiteString(novel_alleles))
    in_what1 <- paste0("germline db ", db_name)
    allele_names1 <- .extract_allele_names(db_fasta_file, in_what1)
    in_what2 <- "'novel_alleles'"
    allele_names2 <- .extract_allele_names(novel_alleles, in_what2)
    clashing_names <- intersect(allele_names1, allele_names2)
    if (length(clashing_names) != 0L) {
        in1string <- paste(clashing_names, collapse=", ")
        stop(wmsg("the following allele names in ", in_what2, " are ",
                  "already present in ", in_what1, ": ", in1string))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### augment_germline_db_[VDJ]()
###

.augment_region_db_success_message <- function(destdir, region_type,
                                               db_fasta_file, novel_fasta)
{
    num_db_alleles <- length(fasta.seqlengths(db_fasta_file))
    num_novel_alleles <- length(fasta.seqlengths(novel_fasta))
    augmented_fasta_file <- get_db_fasta_file(destdir, region_type=region_type)
    total_alleles <- length(fasta.seqlengths(augmented_fasta_file))
    stopifnot(num_db_alleles + num_novel_alleles == total_alleles)

    destdir <- paste0(sub("/*$", "", destdir), "/")
    message("New augmented ", region_type, " germline db successfully ",
            "created in ", destdir)
    message("Number of alleles in augmented db: ",
            num_db_alleles, " + ", num_novel_alleles, " = ", total_alleles)
    message("To use it with igblastn(), do something like:")
    message("")
    message("    igblastn(..., germline_db_", region_type,
            "=\"", destdir, "\")")
    message("")
    message("See '?augment_germline_db_", region_type, "' for ",
            "more information.")
}

### 'novel_alleles' must be the path to a FASTA file (possibly
### gz-compressed), or a named DNAStringSet object.
.augment_region_db <- function(db_name, region_type=VDJ_REGION_TYPES,
                               novel_alleles, destdir=".", overwrite=FALSE,
                               verbose=FALSE)
{
    msg <- c("Some shortcomings were identified in the design of the ",
             "augment_germline_db_[VDJ]() functions so they are now ",
             "deprecated. Please do not use them. They will be replaced ",
             "with a better alternative in future versions of the package.")
    .Deprecated(msg=wmsg(msg))
    check_germline_db_name(db_name)
    db_path <- get_germline_db_path(db_name)
    region_type <- match.arg(region_type)
    db_fasta_file <- get_db_fasta_file(db_path, region_type=region_type)

    novel_fasta <- .novel_alleles_as_temp_edited_fasta(novel_alleles)
    on.exit(unlink(novel_fasta))
    .check_novel_allele_names(db_fasta_file, db_name, novel_fasta)

    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    if (!dir.exists(destdir))
        dir.create(destdir)

    ## Prepare 'fasta_files'. Note that the names we put on 'fasta_files'
    ## will be used by create_region_db() to rename the FASTA files when
    ## they get copied to the <destdir>/<region_type>_original_fasta/ folder.
    fasta_files <- c(db_fasta_file, novel_fasta)
    db_fasta_destfile <- paste0("db--", db_name, "--", region_type, ".fasta")
    names(fasta_files) <- c(db_fasta_destfile, "novel_alleles.fasta")

    ## Create the new region db.
    create_region_db(fasta_files, destdir, region_type=region_type,
                     overwrite=overwrite, verbose=verbose)
    pattern <- paste0("^", region_type, "\\.fasta$")
    make_blastdbs(destdir, pattern=pattern, force=TRUE)

    ## Success!
    .augment_region_db_success_message(destdir, region_type,
                                       db_fasta_file, novel_fasta)
}

### 'novel_alleles' must be the path to a FASTA file or a named
### DNAStringSet object.
augment_germline_db_V <- function(db_name, novel_alleles, destdir=".",
                                  overwrite=FALSE, verbose=FALSE)
{
    .augment_region_db(db_name, region_type="V", novel_alleles,
                       destdir=destdir, overwrite=overwrite, verbose=verbose)
}

augment_germline_db_D <- function(db_name, novel_alleles, destdir=".",
                                  overwrite=FALSE, verbose=FALSE)
{
    .augment_region_db(db_name, region_type="D", novel_alleles,
                       destdir=destdir, overwrite=overwrite, verbose=verbose)
}

augment_germline_db_J <- function(db_name, novel_alleles, destdir=".",
                                  overwrite=FALSE, verbose=FALSE)
{
    .augment_region_db(db_name, region_type="J", novel_alleles,
                       destdir=destdir, overwrite=overwrite, verbose=verbose)
}


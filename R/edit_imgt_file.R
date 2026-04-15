### =========================================================================
### R wrapper to the edit_imgt_file.pl script included in IgBLAST
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_edit_imgt_file_Perl_script()
###

### Requires Perl!
### Checks that Perl script edit_imgt_file.pl is available and that Perl
### is functioning.
get_edit_imgt_file_Perl_script <- function()
{
    igblast_root <- get_igblast_root()
    bin_dir <- get_igblast_root_subdir(igblast_root, "bin")
    script <- file.path(bin_dir, "edit_imgt_file.pl")
    if (!file.exists(script)) {
        details <- c("Perl script 'edit_imgt_file.pl' (needed by ",
                     "internal helper edit_imgt_file()) not found ",
                     "in 'bin' subdirectory.")
        stop_on_invalid_igblast_root(igblast_root, details)
    }
    if (!has_perl())
        stop(wmsg("Setup error: Perl not found."))
    script
}

### BUG: Perl script edit_imgt_file.pl fails to clean the headers of
### most of the FASTA files that IMGT provides for Mus_musculus_C57BL6J.
### As a consequence, installing any of the following germline dbs
### breaks 'list_germline_dbs(long.listing=TRUE)', because most headers
### in the resulting V.fasta/D.fasta/J.fasta files are still the original
### messy 15-field headers:
###   IMGT-202343-3.Mus_musculus_C57BL6J.TRA+TRB
###   IMGT-202405-2.Mus_musculus_C57BL6J.IGH+IGK+IGL
###   IMGT-202405-2.Mus_musculus_C57BL6J.TRA+TRB+TRG+TRD
###   IMGT-202518-3.Mus_musculus_C57BL6J.IGH+IGK+IGL
###   IMGT-202518-3.Mus_musculus_C57BL6J.TRA+TRB+TRG+TRD
### The reason edit_imgt_file.pl fails to clean these headers is because
### they contain only 14 pipes when edit_imgt_file.pl expects exactly 15.
.check_edit_imgt_file_output <- function(outfasta, Perl_script)
{
    allele_names <- names(fasta.seqlengths(outfasta))
    if (any(grepl("|", allele_names, fixed=TRUE))) {
        outfasta <- file_path_as_absolute(outfasta)
        stop(wmsg("script ", Perl_script, " failed to process ",
                  "file ", outfasta, ": some headers still contain ",
                  "the pipe character"))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### edit_imgt_file()
###

### Requires Perl!
edit_imgt_file <- function(infasta, outfasta, errfile=NULL, Perl_script=NULL,
                           check.output=FALSE)
{
    if (!isSingleNonWhiteString(infasta))
        stop(wmsg("'infasta' must be a single (non-empty) string"))
    infasta <- path.expand(infasta)
    if (!isSingleNonWhiteString(outfasta))
        stop(wmsg("'outfasta' must be a single (non-empty) string"))
    outfasta <- path.expand(outfasta)
    if (is.null(errfile)) {
        errfile <- tempfile("edit_imgt_file_errors", fileext=".txt")
        on.exit(unlink(errfile))
    } else if (isSingleNonWhiteString(errfile)) {
        errfile <- path.expand(errfile)
    } else {
        stop(wmsg("'errfile' must be NULL or a single (non-empty) string"))
    }
    if (is.null(Perl_script)) {
        Perl_script <- get_edit_imgt_file_Perl_script()
    } else if (isSingleNonWhiteString(Perl_script)) {
        Perl_script <- path.expand(Perl_script)
    } else {
        stop(wmsg("'Perl_script' must be NULL or a single (non-empty) string"))
    }
    if (!isTRUEorFALSE(check.output))
        stop(wmsg("'check.output' must be TRUE or FALSE"))

    ## This does not work on Windows!
    #system3(Perl_script, outfasta, errfile, args=infasta)

    ## Note that running the Perl script with 'script ...' works on Linux
    ## and Mac but not on Windows. So we run it with 'perl script ...'
    ## instead. This seems to run everywhere.
    system3("perl", outfasta, errfile, args=c(Perl_script, infasta))

    if (check.output)
        .check_edit_imgt_file_output(outfasta, Perl_script)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### redit_imgt_file()
###
### An R implementation of the edit_imgt_file.pl script included in IgBLAST.
###

.edit_imgt_BStringSet_object <- function(dna, in_what)
{
    stopifnot(is(dna, "BStringSet"))
    what <- paste0("some allele names in ", in_what)
    names(dna) <- clean_imgt_fasta_headers(names(dna), what)
    remove_gaps(dna)
}

### Some IMGT FASTA files (e.g. for Aotus_nancymaae and Nonhuman_primates)
### have nucleotide sequences that contain the letter 'x'. Not sure what
### that's supposed to represent. Note that a well established consensus
### is to use 'n' or 'N' to represent an unknown nucleotide (wildcard).
### Anyways, this breaks readDNAStringSet() so we use readBStringSet()
### instead.
redit_imgt_file <- function(infasta, outfasta)
{
    if (!isSingleNonWhiteString(infasta))
        stop(wmsg("'infasta' must be a single (non-empty) string"))
    infasta <- file_path_as_absolute(infasta)
    if (!isSingleNonWhiteString(outfasta))
        stop(wmsg("'outfasta' must be a single (non-empty) string"))
    outfasta <- path.expand(outfasta)

    dna <- readBStringSet(infasta)
    dna <- .edit_imgt_BStringSet_object(dna, infasta)
    writeXStringSet(dna, outfasta)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_redit_imgt_file()
###
### See validate_redit_imgt_file_on_IMGT_release() in
### R/download_IMGT_germline_sequences.R for some extensive validation
### of redit_imgt_file() that we performed on Sep 9, 2025 on various
### IMGT/V-QUEST releases.
###

.compare_edit_imgt_file_vs_redit_imgt_file <- function(fasta_file)
{
    out1 <- tempfile()
    out2 <- tempfile()
    edit_imgt_file(fasta_file, out1)
    redit_imgt_file(fasta_file, out2)
    y1 <- readBStringSet(out1)
    y2 <- readBStringSet(out2)
    identical(as.character(y1), as.character(y2))
}

### Requires Perl!
### Returns number of failures.
validate_redit_imgt_file <- function(dirpath=".", recursive=FALSE)
{
    fasta_files <- list_fasta_files(dirpath, recursive=recursive)
    ## We know that Perl script edit_imgt_file.pl doesn't work properly
    ## on IMGT FASTA files for Mus_musculus_C57BL6J (see BUG above), so
    ## we skip them.
    fasta_files <- grep("/Mus_musculus_C57BL6J/", fasta_files,
                        value=TRUE, invert=TRUE)
    failures <- 0L
    for (i in seq_along(fasta_files)) {
        fasta_file <- fasta_files[[i]]
        message(i, "/", length(fasta_files), " - validate redit_imgt_file() ",
                "on ", fasta_file, " ... ", appendLF=FALSE)
        ok <- .compare_edit_imgt_file_vs_redit_imgt_file(fasta_file)
        if (ok) {
            message("ok")
        } else {
            failures <- failures + 1L
            message("FAILED!")
        }
    }
    failures
}


### =========================================================================
### create_region_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .combine_and_edit_fasta_files()
###

### Workhorse behind create_region_db().
###
### See procedure described at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### for how to create a germline or C-region db from the FASTA files
### available at IMGT. Note that the same procedure can be applied to
### the FASTA files available at AIRR-community/OGRDB.
### This is a 3-step procedure: (1) combine, (2) edit, (3) compile.
### The .combine_and_edit_fasta_files() function below implements
### steps (1) and (2). Perl is required for step (2).
### Compilation (with makeblastdb) will happen at a latter time.
.combine_and_edit_fasta_files <-
    function(fasta_files, destdir, edit_fasta_script,
             region_type=c("V", "D", "J", "C"))
{
    if (!is.character(fasta_files) || anyNA(fasta_files))
        stop(wmsg("'fasta_files' must be a character vector with no NAs"))
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    region_type <- match.arg(region_type)

    ## (1) Combine FASTA files.
    combined_fasta <- file.path(destdir, paste0(".", region_type, ".fasta"))
    concatenate_files(fasta_files, combined_fasta)

    ## (2a) Edit combined FASTA file with 'edit_imgt_file.pl'.
    edited_fasta <- file.path(destdir, paste0(region_type, ".fasta"))
    errfile <- file.path(destdir, paste0(region_type,
                                         "_imgt_script_errors.txt"))

    ## This does not work on Windows!
    #system3(edit_fasta_script, edited_fasta, errfile, args=combined_fasta)

    ## Note that running the Perl script with 'script ...' runs on Linux
    ## and Mac but not on Windows. So we run it with 'perl script ...'
    ## instead. This seems to run everywhere.
    system3("perl", edited_fasta, errfile,
            args=c(edit_fasta_script, combined_fasta))
    unlink(combined_fasta, force=TRUE)

    ## (2b) Mangle seq ids to make them unique if they're not.
    disambiguate_fasta_seqids(edited_fasta)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_region_db()
###

### Perl required!
###
### Creates a "region db" (V-, D-, J-, or C-region) from a collection of
### FASTA files (typically obtained from IMGT or AIRR-community/OGRDB) for
### a given organism.
### See .combine_and_edit_fasta_files() above in this file for the workhorse
### behind create_region_db().
### 'destdir' must be the path to a writable directory that already exists!
### The following subdirectory and files will be added to 'destdir':
###   - V_original_fasta/: subdirectory containing the input FASTA files
###         corresponding to the V regions, one FASTA file per region;
###   - V.fasta: the combined and edited FASTA file produced by calling
###         .combine_and_edit_fasta_files() on the files in V_original_fasta/,
###         with allele names disambiguated if needed.
create_region_db <- function(fasta_files, destdir,
                             region_type=c("V", "D", "J", "C"),
                             edit_fasta_script=NULL)
{
    if (!is.character(fasta_files) || anyNA(fasta_files))
        stop(wmsg("'fasta_files' must be a character vector with no NAs"))
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    region_type <- match.arg(region_type)

    if (is.null(edit_fasta_script)) {
        ## Check that Perl script edit_imgt_file.pl is available and
        ## that Perl is functioning.
        edit_fasta_script <- get_edit_imgt_file_Perl_script()
    }

    ## Create "original fasta" subdir and copy fasta files to it.
    original_fasta_basename <- paste0(region_type, "_original_fasta")
    original_fasta_path <- file.path(destdir, original_fasta_basename)
    stopifnot(dir.create(original_fasta_path))
    stopifnot(all(file.copy(fasta_files, original_fasta_path)))

    ## Combine and edit the original fasta files.
    original_files <- list.files(original_fasta_path, full.names=TRUE)
    .combine_and_edit_fasta_files(original_files, destdir,
                                  edit_fasta_script,
                                  region_type=region_type)
}


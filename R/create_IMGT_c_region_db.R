### =========================================================================
### create_IMGT_c_region_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_and_edit_IMGT_fasta_files()
###

### This is the workhorse behind create_IMGT_c_region_db() and
### create_IMGT_germline_db().
###
### See procedure described at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### for how to create a germline or C-region db from the FASTA files
### available at IMGT.
### This is a 3-step procedure: (1) combine, (2) edit, (3) compile.
### The combine_and_edit_IMGT_fasta_files() function below implements
### steps (1) and (2). Perl is required for step (2).
### Compilation (with makeblastdb) will happen at a latter time.
combine_and_edit_IMGT_fasta_files <-
    function(fasta_files, destdir, edit_fasta_script,
             gene_segment=c("V", "D", "J", "C"))
{
    gene_segment <- match.arg(gene_segment)
    before_edit_dir <- file.path(destdir, "before_edit")
    if (!dir.exists(before_edit_dir))
        dir.create(before_edit_dir, recursive=TRUE)
    unedited_file <- file.path(before_edit_dir, paste0(gene_segment, ".fasta"))
    concatenate_files(fasta_files, unedited_file)
    edited_file <- file.path(destdir, paste0(gene_segment, ".fasta"))
    errfile <- file.path(destdir, paste0(gene_segment,
                                         "_imgt_script_errors.txt"))

    ## This does not work on Windows!
    #system3(edit_fasta_script, edited_file, errfile, args=unedited_file)

    ## Note that running the Perl script with 'script ...' runs on Linux
    ## and Mac but not on Windows. So we run it with 'perl script ...'
    ## instead. This seems to run everywhere.
    system3("perl", edited_file, errfile,
            args=c(edit_fasta_script, unedited_file))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_IMGT_c_region_db()
###

.list_IMGT_C_fasta_files <- function(dirpath)
{
    pattern <- paste0("^IG[HKL]C\\.fasta$")
    files <- list.files(dirpath, pattern=pattern)
    if (length(files) == 0L)
        stop(wmsg("Anomaly: no C-region files found in ", dirpath))
    file.path(dirpath, files)
}

.stop_on_existing_IMGT_c_region_db <- function(destdir)
{
    db_name <- basename(destdir)
    msg1 <- c("C-region db ", db_name, " is already installed.")
    msg2 <- c("Use list_c_region_dbs() to list the C-region databases ",
              "already installed on your machine (see '?list_c_region_dbs').")
    msg3 <- c("Use 'force=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### Perl required!
###
### Create a C-region db from the FASTA files provided by IMGT for a given
### organism. See combine_and_edit_IMGT_fasta_files() above in this file for
### the workhorse behind create_IMGT_c_region_db().
###
### 'organism_path' must be the path to one of the subfolders in
### 'system.file(package="igblastr", "extdata", "constant_regions", "IMGT")'.
###
### 'destdir' will typically be the path to a subdir of
###  <igblastr-cache>/c_region_dbs/
create_IMGT_c_region_db <- function(organism_path, destdir, force=FALSE)
{
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    edit_fasta_script <- get_edit_imgt_file_Perl_script()
    if (dir.exists(destdir) && !force)
        .stop_on_existing_IMGT_c_region_db(destdir)

    ## We first make the db in a temporary folder, and if successful, we
    ## replace 'destdir' with the temporary folder. This achieves atomicity
    ## and avoids loosing the content of the existing 'destdir' in case
    ## something goes wrong.
    tmpdestdir <- tempfile("c_region_db_")
    dir.create(tmpdestdir, recursive=TRUE)
    on.exit(unlink(tmpdestdir, recursive=TRUE, force=TRUE))
    fasta_files <- .list_IMGT_C_fasta_files(organism_path)
    combine_and_edit_IMGT_fasta_files(fasta_files, tmpdestdir,
                                      edit_fasta_script, gene_segment="C")
    replace_file(destdir, tmpdestdir)
}


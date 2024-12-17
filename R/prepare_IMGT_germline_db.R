### =========================================================================
### prepare_IMGT_germline_db()
### -------------------------------------------------------------------------
###
### See procedure described at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### for how to prepare a germline or C-region db from the FASTA files
### available at IMGT.
### This is a 3-step procedure: (1) combine, (2) edit, (3) compile.
### The combine_and_edit_IMGT_fasta_files() function below implements
### steps (1) and (2). Perl is required for step (2).
### Compilation (with makeblastdb) will happen at a latter time.
###
### Nothing in this file is exported.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_and_edit_IMGT_fasta_files()
###
### This is the workhorse behind prepare_IMGT_germline_db() and
### prepare_IMGT_c_region_db().
###

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
### prepare_IMGT_germline_db()
###

.list_IMGT_fasta_files <- function(dirpath, gene_segment=c("V", "D", "J"),
                                   expected_files)
{
    gene_segment <- match.arg(gene_segment)
    pattern <- paste0(gene_segment, "\\.fasta$")
    files <- list.files(dirpath, pattern=pattern)
    if (length(files) == 0L)
        stop(wmsg("Anomaly: no ", gene_segment, " files found in ", dirpath))
    if (!setequal(files, expected_files)) {
        if (all(files %in% expected_files)) {
            warning("some ", gene_segment, " files are missing ",
                    "in ", dirpath, "/ compared to Homo_sapiens")
        } else {
            warning("set of ", gene_segment, " files in ", dirpath, "/ not ",
                    "the same as for Homo_sapiens")
        }
    }
    file.path(dirpath, files)
}

.list_V_files_in_IMGT_IG <- function(IG_path)
{
    EXPECTED_FILES <- paste0("IG", c("HV", "KV", "LV"), ".fasta")
    .list_IMGT_fasta_files(IG_path, "V", EXPECTED_FILES)
}

.list_D_files_in_IMGT_IG <- function(IG_path)
{
    EXPECTED_FILES <- paste0("IGHD", ".fasta")
    .list_IMGT_fasta_files(IG_path, "D", EXPECTED_FILES)
}

.list_J_files_in_IMGT_IG <- function(IG_path)
{
    EXPECTED_FILES <- paste0("IG", c("HJ", "KJ", "LJ"), ".fasta")
    .list_IMGT_fasta_files(IG_path, "J", EXPECTED_FILES)
}

.list_V_files_in_IMGT_TR <- function(TR_path)
{
    EXPECTED_FILES <- paste0("TR", c("AV", "BV", "DV", "GV"), ".fasta")
    .list_IMGT_fasta_files(TR_path, "V", EXPECTED_FILES)
}

.list_D_files_in_IMGT_TR <- function(TR_path)
{
    EXPECTED_FILES <- paste0("TR", c("BD", "DD"), ".fasta")
    .list_IMGT_fasta_files(TR_path, "D", EXPECTED_FILES)
}

.list_J_files_in_IMGT_TR <- function(TR_path)
{
    EXPECTED_FILES <- paste0("TR", c("AJ", "BJ", "DJ", "GJ"), ".fasta")
    .list_IMGT_fasta_files(TR_path, "J", EXPECTED_FILES)
}

.list_V_files_in_IMGT_IG_TR <- function(organism_path)
{
    IG_files <- .list_V_files_in_IMGT_IG(file.path(organism_path, "IG"))
    TR_files <- .list_V_files_in_IMGT_TR(file.path(organism_path, "TR"))
    c(IG_files, TR_files)
}

.list_D_files_in_IMGT_IG_TR <- function(organism_path)
{
    IG_files <- .list_D_files_in_IMGT_IG(file.path(organism_path, "IG"))
    TR_files <- .list_D_files_in_IMGT_TR(file.path(organism_path, "TR"))
    c(IG_files, TR_files)
}

.list_J_files_in_IMGT_IG_TR <- function(organism_path)
{
    IG_files <- .list_J_files_in_IMGT_IG(file.path(organism_path, "IG"))
    TR_files <- .list_J_files_in_IMGT_TR(file.path(organism_path, "TR"))
    c(IG_files, TR_files)
}

.process_IMGT_fasta_files <-
    function(srcdir, destdir, list_files_FUN,
             edit_fasta_script, gene_segment=c("V", "D", "J"))
{
    gene_segment <- match.arg(gene_segment)
    fasta_files <- list_files_FUN(srcdir)
    combine_and_edit_IMGT_fasta_files(fasta_files, destdir,
                                      edit_fasta_script,
                                      gene_segment=gene_segment)
}

.build_IMGT_IG_germline_db <- function(organism_path, destdir,
                                       edit_fasta_script)
{
    IG_path <- file.path(organism_path, "IG")
    .process_IMGT_fasta_files(IG_path, destdir, .list_V_files_in_IMGT_IG,
                              edit_fasta_script, gene_segment="V")
    .process_IMGT_fasta_files(IG_path, destdir, .list_D_files_in_IMGT_IG,
                              edit_fasta_script, gene_segment="D")
    .process_IMGT_fasta_files(IG_path, destdir, .list_J_files_in_IMGT_IG,
                              edit_fasta_script, gene_segment="J")
}

.build_IMGT_TR_germline_db <- function(organism_path, destdir,
                                       edit_fasta_script)
{
    TR_path <- file.path(organism_path, "TR")
    .process_IMGT_fasta_files(TR_path, destdir, .list_V_files_in_IMGT_TR,
                              edit_fasta_script, gene_segment="V")
    .process_IMGT_fasta_files(TR_path, destdir, .list_D_files_in_IMGT_TR,
                              edit_fasta_script, gene_segment="D")
    .process_IMGT_fasta_files(TR_path, destdir, .list_J_files_in_IMGT_TR,
                              edit_fasta_script, gene_segment="J")
}

.build_IMGT_IG_TR_germline_db <- function(organism_path, destdir,
                                          edit_fasta_script)
{
    .process_IMGT_fasta_files(organism_path, destdir,
                              .list_V_files_in_IMGT_IG_TR,
                              edit_fasta_script, gene_segment="V")
    .process_IMGT_fasta_files(organism_path, destdir,
                              .list_D_files_in_IMGT_IG_TR,
                              edit_fasta_script, gene_segment="D")
    .process_IMGT_fasta_files(organism_path, destdir,
                              .list_J_files_in_IMGT_IG_TR,
                              edit_fasta_script, gene_segment="J")
}

.stop_on_existing_IMGT_germline_db <- function(destdir)
{
    db_name <- basename(destdir)
    msg1 <- c("Germline db ", db_name, " is already installed.")
    msg2 <- c("Use list_germline_dbs() to list the germline databases ",
              "already installed on your machine (see '?list_germline_dbs').")
    msg3 <- c("Use 'force=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### Perl required!
###
### 'organism_path' must be the path to an organism subfolder in one of
### the 'IMGT_V-QUEST_reference_directory' folders from the "IMGT-releases"
### store e.g.:
###     <igblastr-cache>
###     └── store
###         └── IMGT-releases
###             └── 202430-2
###                 └── IMGT_V-QUEST_reference_directory
###                     └──  Homo_sapiens
###
### 'destdir' will typically be the path to a subdir of:
###     <igblastr-cache>
###     └── germline_dbs
prepare_IMGT_germline_db <- function(organism_path, destdir,
                                     db_type=c("IG", "TR", "IG-TR"),
                                     force=FALSE)
{
    db_type <- match.arg(db_type)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    edit_fasta_script <- get_edit_imgt_file_Perl_script()
    if (dir.exists(destdir) && !force)
        .stop_on_existing_IMGT_germline_db(destdir)

    ## We check these 2 paths up front so we don't have to check them
    ## later in the .build_IMGT_*_db() functions.
    IG_path <- file.path(organism_path, "IG")
    if (!dir.exists(IG_path))
        stop(wmsg("Anomaly: directory ", IG_path, " not found"))
    TR_path <- file.path(organism_path, "TR")
    if (!dir.exists(TR_path))
        stop(wmsg("Anomaly: directory ", TR_path, " not found"))

    FUN <- switch(db_type,
        "IG"   =.build_IMGT_IG_germline_db,
        "TR"   =.build_IMGT_TR_germline_db,
        "IG-TR"=.build_IMGT_IG_TR_germline_db,
        stop(db_type, ": invalid 'db_type'")
    )

    ## We first make the db in a temporary folder, and if successful, we
    ## replace 'destdir' with the temporary folder. This achieves atomicity
    ## and avoids loosing the content of the existing 'destdir' in case
    ## something goes wrong.
    tmpdestdir <- tempfile("germline_db_")
    dir.create(tmpdestdir, recursive=TRUE)
    on.exit(unlink(tmpdestdir, recursive=TRUE, force=TRUE))
    FUN(organism_path, tmpdestdir, edit_fasta_script)
    replace_file(destdir, tmpdestdir)
}


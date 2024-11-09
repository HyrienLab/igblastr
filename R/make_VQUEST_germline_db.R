### =========================================================================
### make_VQUEST_germline_db()
### -------------------------------------------------------------------------
###
### Execute in R the instructions given at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### to make a V-QUEST germline db. Perl required!
###
### Nothing in this file is exported.


.list_VQUEST_fasta_files <- function(dirpath, group=c("V", "D", "J"),
                                     expected_files)
{
    group <- match.arg(group)
    pattern <- paste0(group, "\\.fasta$")
    files <- list.files(dirpath, pattern=pattern)
    if (length(files) == 0L)
        stop(wmsg("Anomaly: no ", group, " files found in ", dirpath))
    if (!setequal(files, expected_files)) {
        if (all(files %in% expected_files)) {
            warning("some ", group, " files are missing ",
                    "in ", dirpath, "/ compared to Homo_sapiens")
        } else {
            warning("set of ", group, " files in ", dirpath, "/ not ",
                    "the same as for Homo_sapiens")
        }
    }
    file.path(dirpath, files)
}

.list_V_files_in_VQUEST_IG <- function(IG_path)
{
    EXPECTED_FILES <- paste0("IG", c("HV", "KV", "LV"), ".fasta")
    .list_VQUEST_fasta_files(IG_path, "V", EXPECTED_FILES)
}

.list_D_files_in_VQUEST_IG <- function(IG_path)
{
    EXPECTED_FILES <- paste0("IGHD", ".fasta")
    .list_VQUEST_fasta_files(IG_path, "D", EXPECTED_FILES)
}

.list_J_files_in_VQUEST_IG <- function(IG_path)
{
    EXPECTED_FILES <- paste0("IG", c("HJ", "KJ", "LJ"), ".fasta")
    .list_VQUEST_fasta_files(IG_path, "J", EXPECTED_FILES)
}

.list_V_files_in_VQUEST_TR <- function(TR_path)
{
    EXPECTED_FILES <- paste0("TR", c("AV", "BV", "DV", "GV"), ".fasta")
    .list_VQUEST_fasta_files(TR_path, "V", EXPECTED_FILES)
}

.list_D_files_in_VQUEST_TR <- function(TR_path)
{
    EXPECTED_FILES <- paste0("TR", c("BD", "DD"), ".fasta")
    .list_VQUEST_fasta_files(TR_path, "D", EXPECTED_FILES)
}

.list_J_files_in_VQUEST_TR <- function(TR_path)
{
    EXPECTED_FILES <- paste0("TR", c("AJ", "BJ", "DJ", "GJ"), ".fasta")
    .list_VQUEST_fasta_files(TR_path, "J", EXPECTED_FILES)
}

.list_V_files_in_VQUEST_IG_TR <- function(organism_path)
{
    IG_files <- .list_V_files_in_VQUEST_IG(file.path(organism_path, "IG"))
    TR_files <- .list_V_files_in_VQUEST_TR(file.path(organism_path, "TR"))
    c(IG_files, TR_files)
}

.list_D_files_in_VQUEST_IG_TR <- function(organism_path)
{
    IG_files <- .list_D_files_in_VQUEST_IG(file.path(organism_path, "IG"))
    TR_files <- .list_D_files_in_VQUEST_TR(file.path(organism_path, "TR"))
    c(IG_files, TR_files)
}

.list_J_files_in_VQUEST_IG_TR <- function(organism_path)
{
    IG_files <- .list_J_files_in_VQUEST_IG(file.path(organism_path, "IG"))
    TR_files <- .list_J_files_in_VQUEST_TR(file.path(organism_path, "TR"))
    c(IG_files, TR_files)
}

### Note that running the Perl script with 'script ...' runs on Linux and
### Mac but not on Windows. So we run it with 'perl script ...' instead
### which seems to run everywhere.
.edit_VQUEST_fasta <- function(script, infile, outfile, errfile)
{
    ## This does not work on Windows!
    #status <- system2(script, args=infile, stdout=outfile, stderr=errfile)
    status <- system2("perl", args=c(script, infile),
                      stdout=outfile, stderr=errfile)
    errmsg <- readLines(errfile)
    if (length(errmsg) != 0L)
        stop(paste(errmsg, collapse="\n"))
    if (status != 0)
        stop(wmsg("command '", script, " ", infile, "' failed"))
    unlink(errfile)
}

.process_VQUEST_fasta_files <-
    function(srcdir, destdir, list_files_FUN,
             edit_fasta_script, group=c("V", "D", "J"))
{
    group <- match.arg(group)
    files <- list_files_FUN(srcdir)
    before_edit_dir <- file.path(destdir, "before_edit")
    if (!dir.exists(before_edit_dir))
        dir.create(before_edit_dir, recursive=TRUE)
    unedited_file <- file.path(before_edit_dir, paste0(group, ".fasta"))
    concatenate_files(files, unedited_file)
    edited_file <- file.path(destdir, paste0(group, ".fasta"))
    errfile <- file.path(destdir, paste0(group, ".edit_error.txt"))
    .edit_VQUEST_fasta(edit_fasta_script, unedited_file, edited_file, errfile)
}

.build_VQUEST_IG_db <- function(organism_path, db_path, edit_fasta_script)
{
    IG_path <- file.path(organism_path, "IG")
    .process_VQUEST_fasta_files(IG_path, db_path, .list_V_files_in_VQUEST_IG,
                                edit_fasta_script, group="V")
    .process_VQUEST_fasta_files(IG_path, db_path, .list_D_files_in_VQUEST_IG,
                                edit_fasta_script, group="D")
    .process_VQUEST_fasta_files(IG_path, db_path, .list_J_files_in_VQUEST_IG,
                                edit_fasta_script, group="J")
}

.build_VQUEST_TR_db <- function(organism_path, db_path, edit_fasta_script)
{
    TR_path <- file.path(organism_path, "TR")
    .process_VQUEST_fasta_files(TR_path, db_path, .list_V_files_in_VQUEST_TR,
                                edit_fasta_script, group="V")
    .process_VQUEST_fasta_files(TR_path, db_path, .list_D_files_in_VQUEST_TR,
                                edit_fasta_script, group="D")
    .process_VQUEST_fasta_files(TR_path, db_path, .list_J_files_in_VQUEST_TR,
                                edit_fasta_script, group="J")
}

.build_VQUEST_IG_TR_db <- function(organism_path, db_path, edit_fasta_script)
{
    .process_VQUEST_fasta_files(organism_path, db_path,
                                .list_V_files_in_VQUEST_IG_TR,
                                edit_fasta_script, group="V")
    .process_VQUEST_fasta_files(organism_path, db_path,
                                .list_D_files_in_VQUEST_IG_TR,
                                edit_fasta_script, group="D")
    .process_VQUEST_fasta_files(organism_path, db_path,
                                .list_J_files_in_VQUEST_IG_TR,
                                edit_fasta_script, group="J")
}

.stop_on_existing_VQUEST_db <- function(db_name)
{
    msg1 <- c("Germline db ", db_name, " is already installed.")
    msg2 <- c("Use list_germline_dbs() to list the germline databases ",
              "already installed on your machine (see '?list_germline_dbs').")
    msg3 <- c("Use 'force=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

make_VQUEST_germline_db <- function(organism_path, db_path,
                                    db_type=c("IG", "TR", "IG-TR"),
                                    force=FALSE)
{
    db_type <- match.arg(db_type)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    edit_fasta_script <- get_edit_imgt_file_Perl_script()
    if (dir.exists(db_path) && !force)
        .stop_on_existing_VQUEST_db(basename(db_path))

    ## We check these 2 paths up front so we don't have to check them
    ## later in the .build_VQUEST_*_db() functions.
    IG_path <- file.path(organism_path, "IG")
    if (!dir.exists(IG_path))
        stop(wmsg("Anomaly: directory ", IG_path, " not found"))
    TR_path <- file.path(organism_path, "TR")
    if (!dir.exists(TR_path))
        stop(wmsg("Anomaly: directory ", TR_path, " not found"))

    FUN <- switch(db_type,
        "IG"   =.build_VQUEST_IG_db,
        "TR"   =.build_VQUEST_TR_db,
        "IG-TR"=.build_VQUEST_IG_TR_db,
        stop(db_type, ": invalid 'db_type'")
    )
    tmp_db_path <- file.path(dirname(db_path), paste0(".", basename(db_path)))
    dir.create(tmp_db_path, recursive=TRUE)
    on.exit(unlink(tmp_db_path, recursive=TRUE, force=TRUE))
    FUN(organism_path, tmp_db_path, edit_fasta_script)
    replace_file(db_path, tmp_db_path)
}


### =========================================================================
### use_germline_db()
### -------------------------------------------------------------------------


get_germline_dbs_path <- function()
{   
    file.path(R_user_dir("igblastr", "cache"), "germline_dbs")
}

get_germline_db_path <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name), db_name != "USING")
    germline_dbs <- get_germline_dbs_path()
    file.path(germline_dbs, db_name)
}

list_germline_dbs <- function()
{
    germline_dbs <- get_germline_dbs_path()
    if (!dir.exists(germline_dbs))
        return(character(0))
    all_db_names <- setdiff(list.files(germline_dbs), "USING")
    sort(all_db_names)
}

.stop_on_no_installed_germline_db_yet <- function()
{
    msg <- c("You don't have any installed germline database yet. ",
             "Use any of the install_*_germline_db() function (e.g. ",
             "install_VQUEST_germline_db()) to install at least one.")
    stop(wmsg(msg))
}

.get_germline_db <- function()
{
    all_db_names <- list_germline_dbs()
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    germline_dbs <- get_germline_dbs_path()
    using_path <- file.path(germline_dbs, "USING")
    if (!file.exists(using_path)) {
        msg <- c("You haven't selected any germline db to use yet. ",
                 "Please select one with use_germline_db(\"<db_name>\"). ",
                 "See '?use_germline_db' for more information.")
        stop(wmsg(msg))
    }
    db_name <- readLines(using_path)
    if (length(db_name) != 1L)
        stop(wmsg("Anomaly: file '", using_path, "' is corrupted."),
             "\n  ",
             wmsg("File should contain exactly one line. ",
                  "Try to repair with use_germline_db(\"<db_name>\"). ",
                  "See '?use_germline_db' for more information."))
    db_path <- file.path(germline_dbs, db_name)
    if (!dir.exists(db_path))
        stop(wmsg("Anomaly: file '", using_path, "' is invalid."),
             "\n  ",
             wmsg("File content ('", db_name, "') is not the name ",
                  "of an installed germline db. ",
                  "Try to repair with use_germline_db(\"<db_name>\"). ",
                  "See '?use_germline_db' for more information."))
    db_name
}

### Uses the 'makeblastdb' executable provided by NCBI IgBLAST to
### process all FASTA files in the germline db, as instructed at:
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### This produces 10 files for each processed file!
.run_makeblastdb_on_fasta_file <- function(file, makeblastdb_exe)
{
    out_name <- sub("\\.fasta$", "", file)
    args <- c("-parse_seqids", "-dbtype nucl",
              paste("-in", file), paste("-out", out_name))
    out <- suppressWarnings(system2(makeblastdb_exe, args=args,
                                    stdout=TRUE, stderr=TRUE))
    status <- attr(out, "status")
    if (!(is.null(status) || isTRUE(all.equal(status, 0L))))
        stop(wmsg(out))
}

.run_makeblastdb_on_all_fasta_files <- function(db_path)
{
    makeblastdb_exe <- get_igblast_exe("makeblastdb")

    oldwd <- getwd()
    setwd(db_path)
    on.exit(setwd(oldwd))

    fasta_files <- list.files(db_path, pattern="\\.fasta$")
    for (f in fasta_files)
        .run_makeblastdb_on_fasta_file(f, makeblastdb_exe)
}

.stop_on_invalid_db_name <- function(db_name)
{
    msg1 <- c("\"", db_name, "\" is not an installed germline db.")
    msg2 <- c("Use list_germline_dbs() to list the germline databases ",
              "already installed on your machine (see '?list_germline_dbs').")
    msg3 <- c("Note that you can use any of the install_*_germline_db() ",
              "function (e.g. install_VQUEST_germline_db()) to install ",
              "additional germline databases.")
   stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

use_germline_db <- function(db_name=NULL)
{
    if (is.null(db_name))
        return(.get_germline_db())

    ## Check 'db_name'.
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    all_db_names <- list_germline_dbs()
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    if (!(db_name %in% all_db_names))
        .stop_on_invalid_db_name(db_name)

    db_path <- get_germline_db_path(db_name)
    .run_makeblastdb_on_all_fasta_files(db_path)
    using_path <- file.path(germline_dbs, "USING")
    writeLines(db_name, using_path)
    db_name
}


### =========================================================================
### make_blastdbs()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


.infer_region_type_from_fasta_filename <- function(fasta_file)
{
    stopifnot(isSingleNonWhiteString(fasta_file),
              has_suffix(fasta_file, ".fasta"))
    sub("\\.fasta$", "", fasta_file)
}

.make_blastdb_hidden_filename <- function(region_type, suffix)
{
    stopifnot(isSingleNonWhiteString(region_type),
              isSingleNonWhiteString(suffix))
    paste0(".", region_type, "_makeblastdb_", suffix)
}

.make_makeblastdb_version_filename <- function(region_type)
    .make_blastdb_hidden_filename(region_type, "version")

.make_makeblastdb_output_filename <- function(region_type)
    .make_blastdb_hidden_filename(region_type, "output")

.make_makeblastdb_errors_filename <- function(region_type)
    .make_blastdb_hidden_filename(region_type, "errors")

.BLASTDB_SUFFIXES <- c("ndb", "nhr", "nin", "njs", "nog",
                       "nos", "not", "nsq", "ntf", "nto")

.expected_blastdb_filenames <- function(fasta_file)
{
    region_type <- .infer_region_type_from_fasta_filename(fasta_file)
    paste0(region_type, ".", .BLASTDB_SUFFIXES)
}

### A FASTA file is considered to NOT need compilation if the two following
### conditions are satisfied:
###   1. The "makeblastdb version" file associated with the FASTA file
###      is present and contains the same makeblastdb version than
###      reported by makeblastdb_version().
###   2. The expected compilation products are present.
### Note that we don't look at timestamps!
.fasta_file_needs_compilation <- function(db_path, fasta_file)
{
    region_type <- .infer_region_type_from_fasta_filename(fasta_file)
    verfile <- .make_makeblastdb_version_filename(region_type)
    expected_filenames <- c(verfile, .expected_blastdb_filenames(fasta_file))
    paths <- file.path(db_path, expected_filenames)
    if (!all(file.exists(paths)))
        return(TRUE)
    raw_version1 <- readLines(paths[[1L]])
    if (length(raw_version1) == 0L)
        return(TRUE)
    version1 <- raw_version1[[1L]]
    version2 <- makeblastdb_version(raw.version=TRUE)[[1L]]
    !identical(version1, version2)
}

.clean_blastdb_files <- function(db_path, fasta_file)
{
    expected_filenames <- .expected_blastdb_filenames(fasta_file)
    paths <- file.path(db_path, expected_filenames)
    unlink(paths)
}

.check_blastdb_files <- function(db_path, fasta_file)
{
    region_type <- .infer_region_type_from_fasta_filename(fasta_file)
    pattern <- paste0("^", region_type, "\\.n")
    blastdb_files <- list.files(db_path, pattern=pattern)
    if (length(blastdb_files) == 0L)
        stop(wmsg2("no blastdb files found for the ",
                   region_type, "-region db in ", db_path, "/"))
    expected_filenames <- paste0(region_type, ".", .BLASTDB_SUFFIXES)
    if (setequal(blastdb_files, expected_filenames))
        return(TRUE)
    msg1 <- c("Set of blastdb files found in ", db_path, "/ for ",
              "the \"", region_type, "\"-region db is not as expected:")
    expected_in1string <- paste0(expected_filenames, collapse=", ")
    found_in1string <- paste0(blastdb_files, collapse=", ")
    warning(wmsg2(msg1),
            "\n  - expected: ", wmsg2(expected_in1string, margin=14),
            "\n  -    found: ", wmsg2(found_in1string, margin=14))
    FALSE
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .run_makeblastdb_on_fasta_file()
###

### Uses the 'makeblastdb' standalone executable distributed with NCBI IgBLAST
### to "compile" a FASTA file into a db usable with igblastn(). This is the
### last step of the 3-step procedure to create a germline or C-region db from
### a collection of FASTA files. See
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### for more information.
### This "compilation" produces 10 files per FASTA file!
.run_makeblastdb_on_fasta_file <- function(fasta_file, makeblastdb_exe,
                                           verbose=FALSE)
{
    region_type <- .infer_region_type_from_fasta_filename(fasta_file)
    if (verbose)
        message("Making ", region_type, " blast db in ", getwd(), "/ ... ",
                appendLF=FALSE)
    errfile <- .make_makeblastdb_errors_filename(region_type)

    ## Record 'makeblastdb' version in local file.
    verfile <- .make_makeblastdb_version_filename(region_type)
    system3(makeblastdb_exe, verfile, errfile, args="-version")

    ## Make the blast db.
    outfile <- .make_makeblastdb_output_filename(region_type)
    args <- c("-parse_seqids", "-dbtype nucl",
              paste("-in", fasta_file), paste("-out", region_type))
    system3(makeblastdb_exe, outfile, errfile, args=args)

    if (verbose)
        message("ok")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_blastdbs()
###

### Remove the blastdb files produced by make_blastdbs().
clean_blastdbs <- function(db_path, pattern="\\.fasta$")
{
    if (!isSingleNonWhiteString(db_path))
        stop(wmsg("'db_path' must be a single (non-empty) string"))
    if (!dir.exists(db_path))
        stop(wmsg("directory ", db_path, " not found"))

    fasta_files <- list_fasta_files(db_path, pattern, full.names=FALSE)
    for (f in fasta_files)
        .clean_blastdb_files(db_path, f)
    remove_hidden_files(db_path)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_blastdbs()
###

### Returns a named logical vector that indicates the status of all the FASTA
### files in the db. The vector has the file names on it. Status is TRUE if a
### file needs compilation and FALSE otherwise.
.get_fasta_files_statuses <- function(db_path, pattern="\\.fasta$")
{
    fasta_files <- list_fasta_files(db_path, pattern, full.names=FALSE)
    vapply(fasta_files,
        function(f) .fasta_file_needs_compilation(db_path, f),
        logical(1), USE.NAMES=TRUE)
}

### Compiles only the FASTA files that are not already compiled, so it's a
### very fast no-op if all the FASTA files in the db are already compiled.
### Returns the named logical vector obtained with .get_fasta_files_statuses()
### above.
make_blastdbs <- function(db_path, pattern="\\.fasta$",
                          force=FALSE, verbose=FALSE)
{
    if (!isSingleNonWhiteString(db_path))
        stop(wmsg("'db_path' must be a single (non-empty) string"))
    if (!dir.exists(db_path))
        stop(wmsg("directory ", db_path, " not found"))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    statuses <- .get_fasta_files_statuses(db_path, pattern=pattern)
    if (force || any(statuses)) {
        fasta_files <- names(statuses)
        if (!force)
            fasta_files <- fasta_files[statuses]
        makeblastdb_exe <- get_igblast_exe("makeblastdb", check=FALSE)
        db_path <- file_path_as_absolute(db_path)
        oldwd <- getwd()
        setwd(db_path)
        on.exit(setwd(oldwd))
        for (f in fasta_files) {
            if (force || .fasta_file_needs_compilation(db_path, f)) {
                .run_makeblastdb_on_fasta_file(f, makeblastdb_exe,
                                               verbose=verbose)
                .check_blastdb_files(db_path, f)
            }
        }
    }
    statuses
}


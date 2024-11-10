### =========================================================================
### compile_germline_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


.BLASTDB_SUFFIXES <- c("ndb", "nhr", "nin", "njs", "nog",
                       "nos", "not", "nsq", "ntf", "nto")

### TODO: wmsg() was replaced with this in S4Vectors >= 0.45.1 so
### drop .wmsg2() and use wmsg() instead.
.wmsg2 <- function(..., margin=2)
{
    width <- getOption("width") - margin
    paste0(strwrap(paste0(c(...), collapse=""), width=width),
           collapse=paste0("\n", strrep(" ", margin)))
}

.expected_blastdb_filenames <- function(fasta_file)
{
    group <- sub("\\.fasta$", "", fasta_file)
    paste0(group, ".", .BLASTDB_SUFFIXES)
}

### We determine whether a FASTA file needs compilation or not simply
### by looking at the presence of the expected compilation products.
### We don't look at timestamps!
.fasta_file_needs_compilation <- function(db_path, fasta_file)
{
    expected_filenames <- .expected_blastdb_filenames(fasta_file)
    paths <- file.path(db_path, expected_filenames)
    !all(file.exists(paths))
}

.clean_blastdb_files <- function(db_path, fasta_file)
{
    expected_filenames <- .expected_blastdb_filenames(fasta_file)
    paths <- file.path(db_path, expected_filenames)
    unlink(paths)
}

.check_blastdb_files <- function(db_path, fasta_file)
{
    group <- sub("\\.fasta$", "", fasta_file)
    pattern <- paste0("^", group, "\\.n")
    blastdb_files <- list.files(db_path, pattern=pattern)
    if (length(blastdb_files) == 0L)
        stop(.wmsg2("no blastdb files found for \"", group, "\" group ",
                    "in ", db_path, "/"))
    expected_filenames <- paste0(group, ".", .BLASTDB_SUFFIXES)
    if (setequal(blastdb_files, expected_filenames))
        return(TRUE)
    msg1 <- c("Set of blastdb files found in ", db_path, "/ ",
              "for \"", group, "\" group is not as expected:")
    expected_in_1string <- paste0(expected_filenames, collapse=", ")
    found_in_1string <- paste0(blastdb_files, collapse=", ")
    warning(.wmsg2(msg1),
            "\n  - expected: ", .wmsg2(expected_in_1string, margin=14),
            "\n  -    found: ", .wmsg2(found_in_1string, margin=14))
    FALSE
}

### "Compiles" a FASTA file, that is, uses the 'makeblastdb' executable
### provided by NCBI IgBLAST to process all FASTA files in the germline db,
### as instructed at:
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### This produces 10 files per FASTA file!
.run_makeblastdb_on_fasta_file <- function(fasta_file, makeblastdb_exe)
{
    group <- sub("\\.fasta$", "", fasta_file)
    args <- c("-parse_seqids", "-dbtype nucl",
              paste("-in", fasta_file), paste("-out", group))
    out <- suppressWarnings(system2(makeblastdb_exe, args=args,
                                    stdout=TRUE, stderr=TRUE))
    status <- attr(out, "status")
    if (!(is.null(status) || isTRUE(all.equal(status, 0L))))
        stop(wmsg(out))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_germline_db()
###

### Remove the blastdb files produced by compile_germline_db().
clean_germline_db <- function(db_path)
{
    if (!isSingleNonWhiteString(db_path))
        stop(wmsg("'db_path' must be a single (non-empty) string"))
    if (!dir.exists(db_path))
        stop(wmsg("directory ", db_path, " not found"))

    fasta_files <- list.files(db_path, pattern="\\.fasta$")
    for (f in fasta_files)
        .clean_blastdb_files(db_path, f)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compile_germline_db()
###

### Returns a named logical vector that indicates the status of all the FASTA
### files in the db. The vector has the file names on it. Status is TRUE if a
### file needs compilation and FALSE otherwise.
.get_fasta_files_statuses <- function(db_path)
{
    fasta_files <- list.files(db_path, pattern="\\.fasta$")
    statuses <- vapply(fasta_files,
                       function(f) .fasta_file_needs_compilation(db_path, f),
                       logical(1))
    setNames(statuses, fasta_files)
}

### Compiles only the FASTA files that are not already compiled, so it's a
### very fast no-op if all the FASTA files in the db are already compiled.
### Returns the named logical vector obtained with .get_fasta_files_statuses()
### above.
compile_germline_db <- function(db_path)
{
    if (!isSingleNonWhiteString(db_path))
        stop(wmsg("'db_path' must be a single (non-empty) string"))
    if (isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (!dir.exists(db_path))
        stop(wmsg("directory ", db_path, " not found"))

    statuses <- .get_fasta_files_statuses(db_path)
    if (any(statuses)) {
        fasta_files <- names(statuses)[statuses]
        makeblastdb_exe <- get_igblast_exe("makeblastdb")
        oldwd <- getwd()
        setwd(db_path)
        on.exit(setwd(oldwd))
        for (f in fasta_files) {
            if (.fasta_file_needs_compilation(db_path, f)) {
                .run_makeblastdb_on_fasta_file(f, makeblastdb_exe)
                .check_blastdb_files(db_path, f)
            }
        }
    }
    statuses
}


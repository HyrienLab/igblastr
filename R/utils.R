### =========================================================================
### Various general purpose low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


### TODO: Move this to S4Vectors (or BiocBaseUtils).
load_package_gracefully <- function(package, ...)
{
    if (!requireNamespace(package, quietly=TRUE))
        stop("Could not load package ", package, ". Is it installed?\n\n  ",
             wmsg("Note that ", ..., " requires the ", package, " package. ",
                  "Please install it with:"),
             "\n\n    BiocManager::install(\"", package, "\")")
}

### "\xc2\xa0" is some kind of weird white space that sometimes creeps
### in when scrapping dirty HTML documents found on the internet.
is_white_str <- function(x) grepl("^\\s*$", x) | x == "\xc2\xa0"

isSingleNonWhiteString <- function(x) isSingleString(x) && !is_white_str(x)

### Vectorized.
has_prefix <- function(x, prefix)
{
    stopifnot(is.character(x), isSingleString(prefix))
    substr(x, 1L, nchar(prefix)) == prefix
}

### Vectorized.
has_suffix <- function(x, suffix)
{
    stopifnot(is.character(x), isSingleString(suffix))
    x_nc <- nchar(x)
    substr(x, x_nc - nchar(suffix) + 1L, x_nc) == suffix
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level file manipulation
###

### Move 'newfile' to 'oldfile' after nuking 'oldfile' if needed.
### Works with files or directories.
replace_file <- function(oldfile, newfile)
{
    stopifnot(isSingleNonWhiteString(oldfile), isSingleNonWhiteString(newfile))
    if (!file.exists(newfile))
        stop(wmsg(newfile, ": no such file or directory"))
    unlink(oldfile, recursive=TRUE, force=TRUE)
    ok <- file.rename(newfile, oldfile)
    if (!ok)
        stop(wmsg("failed to replace '", oldfile, "' with '", newfile, "'"))
}

### Does not recursively search hidden files i.e. it only removes the
### hidden files (and directories if 'include.hidden.dirs=TRUE') that
### are located **directly** under 'path'.
remove_hidden_files <- function(path=".", include.hidden.dirs=FALSE)
{
    stopifnot(isSingleNonWhiteString(path), dir.exists(path),
              isTRUEorFALSE(include.hidden.dirs))
    hidden_files <- list.files(path, pattern="^\\.", all.files=TRUE,
                               full.names=TRUE, no..=TRUE)
    unlink(hidden_files, recursive=include.hidden.dirs, force=TRUE)
}

concatenate_files <- function(files, out=stdout(), n=50000L)
{
    stopifnot(is.character(files))
    if (is.character(out)) {
        out <- file(out, "wb")
        on.exit(close(out))
    }
    for (f in files) {
        con <- file(f, "rb")
        while (TRUE) {
            bytes <- readBin(con, what=raw(), n=n)
            if (length(bytes) == 0L)
                break
            writeBin(bytes, out)
        }
        close(con)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastr_cache()
###

igblastr_cache <- function()
{
    getOption("igblastr_cache", R_user_dir("igblastr", "cache"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### disambiguate_fasta_seqids()
###

### Similar to base::make.unique() but mangles with suffixes made of
### lowercase letters.
.make_pool_of_suffixes <- function(min_pool_size)
{
    max_pool_size <- (length(letters)**8 - 1) / (length(letters) - 1) - 1
    if (min_pool_size > max_pool_size)
        stop(wmsg("too many duplicate seq ids"))
    ans <- character(0)
    for (i in 1:7) {
        ans <- c(ans, mkAllStrings(letters, i))
        if (length(ans) >= min_pool_size)
            return(ans)
    }
    ## Should never happen because we checked for this condition earlier (see
    ## above).
    stop(wmsg("too many duplicate seq ids"))
}

.make_unique_seqids <- function(seqids)
{
    stopifnot(is.character(seqids))
    if (length(seqids) <= 1L)
        return(seqids)
    oo <- order(seqids)
    seqids2 <- seqids[oo]
    ir <- IRanges(1L, runLength(Rle(seqids2)))
    pool_of_suffixes <- .make_pool_of_suffixes(max(width(ir)))
    suffixes <- extractList(pool_of_suffixes, ir)  # CharacterList
    suffixes[lengths(suffixes) == 1L] <- ""
    unlist(suffixes, use.names=FALSE)
    seqids2 <- paste0(seqids2, unlist(suffixes, use.names=FALSE))
    ans <- seqids2[S4Vectors:::reverseIntegerInjection(oo, length(oo))]
    setNames(ans, names(seqids))
}

### In-place replacement!
disambiguate_fasta_seqids <- function(filepath)
{
    stopifnot(isSingleNonWhiteString(filepath))
    fasta_lines <- readLines(filepath)
    header_idx <- grep("^>", fasta_lines)
    header_lines <- fasta_lines[header_idx]
    if (anyDuplicated(header_lines)) {
        fasta_lines[header_idx] <- .make_unique_seqids(header_lines)
        writeLines(fasta_lines, filepath)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### tabulate_c_region_seqids_by_locus()
### tabulate_germline_seqids_by_group()
###

GENE_LOCI <- c("IGH", "IGK", "IGL")

### Returns a named integer vector with GENE_LOCI as names.
tabulate_c_region_seqids_by_locus <- function(seqids)
{
    stopifnot(is.character(seqids))
    prefixes <- substr(seqids, 1L, 3L)
    m <- match(prefixes, GENE_LOCI)
    if (anyNA(m)) {
        in1string <- paste0(GENE_LOCI, collapse=", ")
        stop(wmsg("not all seq ids start with one of the following prefixes: ",
                  in1string))
    }
    setNames(tabulate(m, length(GENE_LOCI)), GENE_LOCI)
}

### Group names are formed by concatenating a locus name (IGH, IGK, or IGL)
### and a region type (V, D, or J).
GERMLINE_GROUPS <- c(paste0("IGH", c("V", "D", "J")),
                     paste0("IGK", c("V", "J")),
                     paste0("IGL", c("V", "J")))

### Returns a named integer vector with GERMLINE_GROUPS as names.
tabulate_germline_seqids_by_group <- function(seqids)
{
    stopifnot(is.character(seqids))
    prefixes <- substr(seqids, 1L, 4L)
    m <- match(prefixes, GERMLINE_GROUPS)
    if (anyNA(m)) {
        in1string <- paste0(GERMLINE_GROUPS, collapse=", ")
        stop(wmsg("not all seq ids start with one of the following prefixes: ",
                  in1string))
    }
    setNames(tabulate(m, length(GERMLINE_GROUPS)), GERMLINE_GROUPS)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_version_file()
###

read_version_file <- function(dirpath)
{
    stopifnot(isSingleNonWhiteString(dirpath))
    version_path <- file.path(dirpath, "version")
    if (!file.exists(version_path))
        stop(wmsg("missing 'version' file in ", dirpath, "/"))
    version <- readLines(version_path)
    if (length(version) != 1L)
        stop(wmsg("file '", version_path, "' should contain exactly one line"))
    version <- trimws(version)
    if (version == "")
        stop(wmsg("file '", version_path, "' contains only white spaces"))
    version
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_db_in_use()
###

### Returns "" if no db is currently in use.
get_db_in_use <- function(dbs_path, what=c("germline", "C-region"))
{
    stopifnot(isSingleNonWhiteString(dbs_path), dir.exists(dbs_path))
    what <- match.arg(what)
    if (what == "germline") {
        fun <- "use_germline_db"
    } else {
        fun <- "use_c_region_db"
    }
    repair_with <- paste0("Try to repair with ", fun, "(\"<db_name>\").")
    see <- paste0("See '?", fun, "' for more information.")

    using_path <- file.path(dbs_path, "USING")
    if (!file.exists(using_path))
        return("")
    db_name <- readLines(using_path)
    if (length(db_name) != 1L)
        stop(wmsg("Anomaly: file '", using_path, "' is corrupted."),
             "\n  ",
             wmsg("File should contain exactly one line. ",
                  repair_with, " ", see))
    db_path <- file.path(dbs_path, db_name)
    if (!dir.exists(db_path))
        stop(wmsg("Anomaly: file '", using_path, "' is invalid."),
             "\n  ",
             wmsg("File content ('", db_name, "') is not the name ",
                  "of a cached ", what, " db. ",
                  repair_with, " ", see))
    db_path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous stuff
###

urlExists <- function(url)
{
    response <- try(HEAD(url, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response), "  Please check your internet connection.")
    response$status_code != 404L
}

getUrlContent <- function(url, query=list(), type=NULL, encoding=NULL)
{
    stopifnot(is.list(query))
    if (length(query) != 0L)
        stopifnot(!is.null(names(query)))
    response <- try(GET(url, user_agent("igblastr"), query=query), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response), "  Please check your internet connection.")
    if (response$status_code == 404L)
        stop(wmsg("Not Found (HTTP 404): ", url))
    stop_for_status(response)
    content(response, type=type, encoding=encoding)
}

### Note that in case of user interrupt (CTRL+C) we can end up with a
### partial download and corrupted file! Can we achieve atomic behavior?
### TODO: Try to make the behavior atomic, under any circumstance.
download_as_tempfile <- function(dir_url, filename, ...)
{
    url <- paste0(dir_url, filename)
    destfile <- tempfile()
    code <- try(suppressWarnings(download.file(url, destfile, ...)),
                silent=TRUE)
    if (inherits(code, "try-error") || code != 0L)
        stop(wmsg("Failed to download ", filename, " ",
                  "from ", dir_url, ". ",
                  "Are you connected to the internet?"))
    destfile
}

### A thin wrapper to untar() with more user-friendly error handling.
### 'exdir' should be the path to an existing directory that is
### preferrably empty.
untar2 <- function(tarfile, original_tarball_name, exdir=".")
{
    stopifnot(isSingleNonWhiteString(tarfile),
              isSingleNonWhiteString(original_tarball_name),
              isSingleNonWhiteString(exdir),
              dir.exists(exdir))
    code <- suppressWarnings(untar(tarfile, exdir=exdir))
    if (code != 0L)
        stop(wmsg("Anomaly: something went wrong during ",
                  "extraction of '", tarfile, "' (the local copy ",
                  "of '", original_tarball_name, "') to '", exdir, "'."))
}

### Returns the OS (e.g. Linux, Windows, or Darwin) and arch (e.g. x86_64
### or arm64) in a character vector of length 2, with names "OS" and "arch".
### Note that if the OS or arch cannot be obtained with Sys.info() then they
### get replaced with an NA.
get_OS_arch <- function()
{
    sys_info <- Sys.info()
    sysname <- sys_info[["sysname"]]
    if (!isSingleNonWhiteString(sysname))
        sysname <- NA_character_
    machine <- sys_info[["machine"]]
    if (!isSingleNonWhiteString(machine))
        machine <- NA_character_
    c(OS=sysname, arch=machine)
}

add_exe_suffix_on_Windows <- function(files, OS=get_OS_arch()[["OS"]])
{
    stopifnot(is.character(files), isSingleStringOrNA(OS))
    if (length(files) == 0L || is.na(OS) || !grepl("^win", tolower(OS)))
        return(files)
    paste0(files, ".exe")
}

system_command_works <- function(command, args=character())
{
    out <- try(suppressWarnings(system2(command, args=args,
                                        stdout=TRUE, stderr=TRUE)),
               silent=TRUE)
    if (inherits(out, "try-error"))
        return(FALSE)
    status <- attr(out, "status")
    is.null(status) || isTRUE(all.equal(status, 0L))
}

has_perl <- function() system_command_works("perl", args="-v")

system3 <- function(command, outfile, errfile, args=character())
{
    status <- system2(command, args=args, stdout=outfile, stderr=errfile)
    errmsg <- readLines(errfile)
    if (length(errmsg) != 0L)
        stop(paste(errmsg, collapse="\n"))
    unlink(errfile)
    if (status != 0) {
        cmd_in_1string <- paste(c(command, args), collapse=" ")
        stop(wmsg("command '", cmd_in_1string, "' failed"))
    }
}

display_local_file_in_browser <- function(file)
{
    top_html <- tempfile()
    writeLines("<PRE>", top_html)
    bottom_html <- tempfile()
    writeLines("</PRE>", bottom_html)
    temp_html <- tempfile(fileext=".html")
    concatenate_files(c(top_html, file, bottom_html), out=temp_html)
    temp_url <- paste0("file://", temp_html)
    browseURL(temp_url)
}


### =========================================================================
### Various general purpose low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


### "\xc2\xa0" is some kind of weird white space that sometimes creeps
### in when scrapping dirty HTML documents found on the internet.
is_white_str <- function(x) grepl("^\\s*$", x) | x == "\xc2\xa0"

isSingleNonWhiteString <- function(x) isSingleString(x) && !is_white_str(x)

### Vectorized.
has_suffix <- function(x, suffix)
{
    stopifnot(is.character(x), isSingleString(suffix))
    x_nc <- nchar(x)
    substr(x, x_nc - nchar(suffix) + 1L, x_nc) == suffix
}

urlExists <- function(url)
{
    response <- try(HEAD(url, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response), "  Please check your internet connection.")
    response$status_code != 404L
}

getUrlContent <- function(url, type=NULL, encoding=NULL)
{
    response <- try(GET(url, user_agent("igblastr")), silent=TRUE)
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


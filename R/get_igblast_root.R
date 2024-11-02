### =========================================================================
### get_igblast_root()
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers used by set_internal_igblast_root(),
### get_igblast_root(), and get_igblast_exe()
### 

.stop_on_invalid_installation <- function(details, igblast_root)
{
    msg <- c("Invalid IgBlast installation at '", igblast_root, "'")
    obtained_via <- attr(igblast_root, "obtained_via")
    if (is.null(obtained_via)) {
        msg <- c("Anomaly: ", msg)
    } else {
        msg <- c("Setup error: ", msg, " ",
                 "(path obtained with '", obtained_via, "')")

    }
    stop(wmsg(msg), ". ", wmsg(details))
}

### Returns "<igblast_root>/bin". Guaranteed to return a valid path.
.get_igblast_root_bin <- function(igblast_root)
{
    if (!dir.exists(igblast_root))
        .stop_on_invalid_installation("Directory does not exist.",
                                      igblast_root)
    bin_dir <- file.path(igblast_root, "bin")
    if (!dir.exists(bin_dir))
        .stop_on_invalid_installation("Directory has no 'bin' subdirectory.",
                                      igblast_root)
    bin_dir
}

### Guaranteed to return the path to a working executable.
.make_igblast_cmd_path <- function(igblast_root, cmd=c("igblastn", "igblastp"))
{
    bin_dir <- .get_igblast_root_bin(igblast_root)
    cmd <- match.arg(cmd)
    cmd_path <- file.path(bin_dir, cmd)
    if (!file.exists(cmd_path)) {
        details <- c("No '", cmd, "' command in 'bin' subdirectory.")
        .stop_on_invalid_installation(details, igblast_root)
    }
    if (!system_command_works(cmd_path, "-version")) {
        details <- c("'", cmd_path, " -version' does not work.")
        .stop_on_invalid_installation(details, igblast_root)
    }
    cmd_path
}

.check_igblast_installation <- function(igblast_root)
{
    bin_dir <- .get_igblast_root_bin(igblast_root)
    ## Check content of 'bin_dir'.
    bin_files <- list.files(bin_dir)
    required_bin_files <- c("igblastn", "igblastp",
                            "makeblastdb", "edit_imgt_file.pl")
    for (file in required_bin_files) {
        if (!(file %in% bin_files)) {
            details <- c("No '", file, "' file in 'bin' subdirectory.")
            .stop_on_invalid_installation(details, igblast_root)
        }
    }
    ## We ignore the returned path. Only purpose is to check that the
    ## igblastn and igblastp executables work.
    .make_igblast_cmd_path(igblast_root, cmd="igblastn")
    .make_igblast_cmd_path(igblast_root, cmd="igblastn")
}

.what_to_do_if_igblast_not_found <- function()
{
    paste0("Please use install_igblast() to download and ",
           "install a pre-compiled IgBlast from NCBI FTP site. ",
           "See '?install_igblast' for the details. ",
           "Alternatively you can set environment variable ",
           "IGBLAST_ROOT to point to an existing IgBlast ",
           "installation on your machine. See '?IGBLAST_ROOT' ",
           "for more information.")
}

### Returns the path to an "external" (i.e. not igblastr-controlled)
### IgBlast installation, if any. Note that the path is returned with
### the "obtained_via" attribute on it. The attribute is a single string
### describing where the path is coming from, that is, whether it's
### coming from global option 'igblast_root' or from environment
### variable IGBLAST_ROOT.
### Checks the returned installation.
.get_external_igblast_root <- function()
{
    ## Looking for an external installation.
    igblast_root <- getOption("igblast_root")
    if (is.null(igblast_root)) {
        igblast_root <- Sys.getenv("IGBLAST_ROOT")
        if (!nzchar(igblast_root))
            stop("No IgBlast installation found.\n  ",
                 wmsg(.what_to_do_if_igblast_not_found()))
        obtained_via <- "Sys.getenv(\"IGBLAST_ROOT\")"
        if (!isSingleNonWhiteString(igblast_root))
            stop(wmsg("Anomaly: '", obtained_via, "' should ",
                      "be a single (non-empty) string"))
    } else {
        obtained_via <- "getOption(\"igblast_root\")"
        if (!isSingleNonWhiteString(igblast_root))
            stop(wmsg("Anomaly: '", obtained_via, "' should ",
                      "be NULL or a single (non-empty) string"))
    }
    attr(igblast_root, "obtained_via") <- obtained_via
    .check_igblast_installation(igblast_root)
    igblast_root  # path has "obtained_via" attribute
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set or get the "internal" (i.e. igblastr-controlled) IgBlast
### installation to use
###

igblastr_local_executables_dir <- function()
{
    file.path(R_user_dir("igblastr", "cache"), "executables")
}

### Sets the path to the "internal" IgBlast installation to use and
### returns it. Only used by install_igblast() at the moment.
### Checks the returned installation.
set_internal_igblast_root <- function(rootbasename)
{
    if (!isSingleNonWhiteString(rootbasename))
        stop(wmsg("'rootbasename' must be a single (non-empty) string"))
    local_executables_dir <- igblastr_local_executables_dir()
    if (!dir.exists(local_executables_dir))
        stop(wmsg("Anomaly: no '", local_executables_dir, "'."))
    igblast_root <- file.path(local_executables_dir, rootbasename)
    .check_igblast_installation(igblast_root)
    using_path <- file.path(local_executables_dir, "USING")
    writeLines(rootbasename, using_path)
    igblast_root
}

### Exported!
use_igblast <- set_internal_igblast_root

### Returns the path to the "internal" IgBlast installation currently
### in use, if any. Otherwise, returns NA_character_.
### NOTE: Unlike set_internal_igblast_root() and .get_external_igblast_root()
### above, this function does NOT check the returned installation.
get_internal_igblast_root <- function()
{
    local_executables_dir <- igblastr_local_executables_dir()
    #if (!dir.exists(local_executables_dir))
    #    return(NA_character_)
    using_path <- file.path(local_executables_dir, "USING")
    if (!file.exists(using_path))
        return(NA_character_)
    rootbasename <- readLines(using_path)
    file.path(local_executables_dir, rootbasename)  # naked path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_igblast_root() and get_igblast_exe()
###

### Returns the path to the IgBlast installation used by igblastr. In case
### of an external installation, the path is returned with the "obtained_via"
### attribute on it. See .get_external_igblast_root() above in this file for
### more information.
### Exported!
get_igblast_root <- function()
{
    igblast_root <- get_internal_igblast_root()
    if (!is.na(igblast_root))
        return(igblast_root)      # naked path
    .get_external_igblast_root()  # path has "obtained_via" attribute
}

get_igblast_exe <- function(cmd=c("igblastn", "igblastp"))
{
    igblast_root <- get_igblast_root()
    cmd <- match.arg(cmd)
    .make_igblast_cmd_path(igblast_root, cmd=cmd)
}


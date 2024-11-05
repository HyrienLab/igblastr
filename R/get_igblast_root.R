### =========================================================================
### get_igblast_root() and set_igblast_root()
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

.required_igblast_root_bin_files <- function(OS)
{
    files <- c("igblastn", "igblastp", "makeblastdb")
    files <- add_exe_suffix_on_Windows(files, OS=OS)
    c(files, "edit_imgt_file.pl")
}

### Guaranteed to return the path to a working executable.
.make_igblast_cmd_path <- function(igblast_root,
                                   cmd=c("igblastn", "igblastp"),
                                   OS=get_OS_arch()[["OS"]])
{
    bin_dir <- .get_igblast_root_bin(igblast_root)
    cmd <- add_exe_suffix_on_Windows(match.arg(cmd), OS=OS)
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

### Returns 'igblast_root' if IgBlast installation is valid. Otherwise raises
### an error.
.check_igblast_installation <- function(igblast_root)
{
    bin_dir <- .get_igblast_root_bin(igblast_root)
    ## Check content of 'bin_dir'.
    bin_files <- list.files(bin_dir)
    OS <- get_OS_arch()[["OS"]]
    required_bin_files <- .required_igblast_root_bin_files(OS)
    for (file in required_bin_files) {
        if (!(file %in% bin_files)) {
            details <- c("No '", file, "' file in 'bin' subdirectory.")
            .stop_on_invalid_installation(details, igblast_root)
        }
    }
    ## We ignore the returned path. Only purpose is to check that the
    ## igblastn and igblastp executables work.
    .make_igblast_cmd_path(igblast_root, cmd="igblastn", OS=OS)
    .make_igblast_cmd_path(igblast_root, cmd="igblastp", OS=OS)
    igblast_root
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set or get the "external" (i.e. not igblastr-managed) IgBlast
### installation to use
###

### Returns the path to an "external" IgBlast installation, if any.
### Otherwise, returns NA_character_.
### Note that the path is returned with the "obtained_via" attribute on it.
### The value of the attribute is a single string describing how the path
### was obtained, that is, whether it's coming from global
### option 'igblast_root' or from environment variable IGBLAST_ROOT.
### Checks the returned installation.
.get_external_igblast_root <- function()
{
    igblast_root <- getOption("igblast_root")
    if (is.null(igblast_root)) {
        igblast_root <- Sys.getenv("IGBLAST_ROOT")
        if (!nzchar(igblast_root))
            return(NA_character_)
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
}

### Checks the returned installation.
.set_external_igblast_root <- function(path)
{
    igblast_root <- file_path_as_absolute(path)
    .check_igblast_installation(igblast_root)
    options(igblast_root=igblast_root)
    igblast_root
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set or get the "internal" (i.e. igblastr-managed) IgBlast
### installation to use
###

igblastr_local_executables_dir <- function()
{
    file.path(R_user_dir("igblastr", "cache"), "executables")
}

### Sets the path to the "internal" IgBlast installation to use and
### returns it.
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
### get_igblast_root() and set_igblast_root()
###

### List all "internal" IgBlast installations ordered by decreasing version.
.list_igblast_rootbasenames <- function()
{
    local_executables_dir <- igblastr_local_executables_dir()
    if (!dir.exists(local_executables_dir))
        return(character(0))
    all_rootbasenames <- setdiff(list.files(local_executables_dir), "USING")
    pattern <- "^ncbi-igblast-([.0-9]+)"
    versions <- numeric_version(sub(pattern, "\\1", all_rootbasenames))
    all_rootbasenames[order(versions, decreasing=TRUE)]
}

### Returns the path to the IgBlast installation used by igblastr. In case
### of an external installation, the path is returned with the "obtained_via"
### attribute on it. See .get_external_igblast_root() above in this file for
### more information.
### Checks the returned installation only if it's an "external" one.
### Exported!
get_igblast_root <- function()
{
    ## 1. Look for an "internal" IgBlast installation that is in use, and
    ##    return it if any.
    igblast_root <- get_internal_igblast_root()
    if (!is.na(igblast_root))
        return(igblast_root)  # naked path

    ## 2. Look for an "external" IgBlast installation and return it if any.
    igblast_root <- .get_external_igblast_root()
    if (!is.na(igblast_root))
        return(igblast_root)  # path has "obtained_via" attribute

    ## 3. Look for "internal" IgBlast installations, choose the first one (i.e.
    ##    the one with the highest version) if there are any, set it as the one
    ##    in use, and return it.
    all_rootbasenames <- .list_igblast_rootbasenames()
    if (length(all_rootbasenames) != 0L) {
        rootbasename <- all_rootbasenames[[1L]]
        return(set_internal_igblast_root(rootbasename))
    }

    ## All the above have failed!
    stop("No IgBlast installation found.\n  ",
         wmsg("Please use install_igblast() to download and ",
              "install a pre-compiled IgBlast from NCBI FTP site. ",
              "See '?install_igblast' for the details. ",
              "Alternatively you can set environment variable ",
              "IGBLAST_ROOT to point to an existing IgBlast ",
              "installation on your machine. See '?IGBLAST_ROOT' ",
              "for more information."))
}

### Sets the path to the IgBlast installation to use and returns it.
### This can be either an "internal" or an "external" installation.
### In the former case, 'path' must be the "root basename" of an
### existing internal installation, that is, its path relative to
### igblastr_local_executables_dir().
### In the latter case, it must be the full path (absolute or relative)
### to the root directory of an existing external installation.
### Always checks the returned installation.
### Exported!
set_igblast_root <- function(path)
{
    if (!isSingleNonWhiteString(path))
        stop(wmsg("'path' must be a single (non-empty) string"))
    if (dir.exists(path)) {
        ## Set the path to the "external" IgBlast installation to use.
        igblast_root <- .set_external_igblast_root(path)
        ## Remove the USING file if any.
        local_executables_dir <- igblastr_local_executables_dir()
        if (dir.exists(local_executables_dir)) {
            using_path <- file.path(local_executables_dir, "USING")
            if (file.exists(using_path))
                unlink(using_path)
        }
        return(igblast_root)
    }
    ## Set the path to the "internal" IgBlast installation to use.
    all_rootbasenames <- .list_igblast_rootbasenames()
    if (length(all_rootbasenames) == 0L) {
        err_msg2 <- c("If there's no IgBlast installation on your ",
                      "machine that you can point 'path' at, then use ",
                      "install_igblast() to download and install a ",
                      "pre-compiled IgBlast from NCBI FTP site. ",
                      "See '?install_igblast' for the details.")
    } else {
        if (path %in% all_rootbasenames)
            return(set_internal_igblast_root(path))
        all_in_1string <- paste0("\"", all_rootbasenames, "\"", collapse=", ")
        err_msg2 <- c("Valid \"root basenames\": ", all_in_1string, ".")
    }
    err_msg1 <- c("'path' must be either the path to an existing directory ",
                  "that contains a valid externally managed IgBlast ",
                  "installation, or it must be a valid \"root basename\" ",
                  "i.e. the \"root basename\" of an existing igblastr-managed ",
                  "IgBlast installation.")
    stop(wmsg(err_msg1), "\n  ", wmsg(err_msg2))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_igblast_exe()
###

get_igblast_exe <- function(cmd=c("igblastn", "igblastp"))
{
    igblast_root <- get_igblast_root()
    cmd <- match.arg(cmd)
    .make_igblast_cmd_path(igblast_root, cmd=cmd)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblast_version() and igblast_info()
###

.extract_version_from_cmd_output <- function(output)
{
    sub("^igblastn: *", "", output[[1L]])
}

### Exported!
igblast_version <- function()
{
    igblastn_exe <- get_igblast_exe("igblastn")
    out <- system2(igblastn_exe, "-version", stdout=TRUE)
    .extract_version_from_cmd_output(out)
}

### Exported!
print.igblast_info <- function(x, ...)
{
    x <- lapply(x, function(x) paste(x, collapse="; "))
    x <- paste0(names(x), ": ", as.character(x))
    cat(x, sep="\n")
}

### Exported!
igblast_info <- function()
{
    igblast_root <- get_igblast_root()
    OS_arch <- get_OS_arch()
    OS <- OS_arch[["OS"]]
    igblastn_exe <- .make_igblast_cmd_path(igblast_root, cmd="igblastn", OS=OS)
    igblastn_version <- system2(igblastn_exe, "-version", stdout=TRUE)
    version <- .extract_version_from_cmd_output(igblastn_version)
    #igblastp_exe <- .make_igblast_cmd_path(igblast_root, cmd="igblastp", OS=OS)
    #igblastp_version <- system2(igblastp_exe, "-version", stdout=TRUE)
    build <- igblastn_version[[2L]]  # should be same as igblastp_version[[2L]]
    build <- sub("^ *Package: ", "", build)

    ans <- list(
        igblast_root=igblast_root,
        Version=version,
        `OS/arch`=paste(OS_arch, collapse="/"),
        Build=build
        #igblastn_exe=igblastn_exe,
        #igblastn_version=igblastn_version[[1L]],
        #igblastp_exe=igblastp_exe,
        #igblastp_version=igblastp_version[[1L]]
    )
    class(ans) <- "igblast_info"
    ans
}


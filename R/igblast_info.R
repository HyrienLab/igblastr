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
    msg <- c("Invalid IgBLAST installation at '", igblast_root, "'")
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

### Returns 'igblast_root' if IgBLAST installation is valid. Otherwise raises
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
### Set or get the "external" (i.e. not igblastr-managed) IgBLAST
### installation to use
###

### Returns the path to an "external" IgBLAST installation, if any.
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
### Set or get the "internal" (i.e. igblastr-managed) IgBLAST
### installation to use
###

get_internal_igblast_roots <- function()
{
    file.path(R_user_dir("igblastr", "cache"), "igblast_roots")
}

### Sets the path to the "internal" IgBLAST installation to use and
### returns it.
### Checks the returned installation.
set_internal_igblast_root <- function(version)
{
    if (!isSingleNonWhiteString(version))
        stop(wmsg("'version' must be a single (non-empty) string"))
    internal_roots <- get_internal_igblast_roots()
    if (!dir.exists(internal_roots))
        stop(wmsg("Anomaly: no '", internal_roots, "'."))
    igblast_root <- file.path(internal_roots, version)
    .check_igblast_installation(igblast_root)
    using_path <- file.path(internal_roots, "USING")
    writeLines(version, using_path)
    igblast_root
}

### Returns the path to the "internal" IgBLAST installation currently
### in use, if any. Otherwise, returns NA_character_.
### NOTE: Unlike set_internal_igblast_root() and .get_external_igblast_root()
### above, this function does NOT check the returned installation.
get_internal_igblast_root <- function()
{
    internal_roots <- get_internal_igblast_roots()
    #if (!dir.exists(internal_roots))
    #    return(NA_character_)
    using_path <- file.path(internal_roots, "USING")
    if (!file.exists(using_path))
        return(NA_character_)
    version <- readLines(using_path)
    if (length(version) != 1L)
        stop(wmsg("Anomaly: '", using_path, "' is corrupted."),
             "\n  ",
             wmsg("File should contain exactly one line. ",
                  "Try to repair with set_igblast_root() ",
                  "(see '?set_igblast_root' for more information)."))
    igblast_root <- file.path(internal_roots, version)
    if (!dir.exists(igblast_root))
        stop(wmsg("Anomaly: '", using_path, "' is invalid."),
             "\n  ",
             wmsg("File content ('", version, "') is not the version ",
                  "of an igblastr-managed installation of IgBLAST. ",
                  "Try to repair with set_igblast_root() ",
                  "(see '?set_igblast_root' for more information)."))
    igblast_root  # naked path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_igblast_root() and set_igblast_root()
###

### List all "internal" IgBLAST installations ordered by decreasing version.
.list_installed_igblast_versions <- function()
{
    internal_roots <- get_internal_igblast_roots()
    if (!dir.exists(internal_roots))
        return(character(0))
    all_versions <- setdiff(list.files(internal_roots), "USING")
    oo <- order(numeric_version(all_versions, strict=FALSE), decreasing=TRUE)
    all_versions[oo]
}

.what_to_do_if_no_internal_installation_yet <- function()
    c("Please use install_igblast() to download and ",
      "install a pre-compiled IgBLAST from NCBI FTP site. ",
      "See '?install_igblast' for the details. ",
      "Alternatively you can set environment variable ",
      "IGBLAST_ROOT to point to an existing IgBLAST ",
      "installation on your machine. See '?IGBLAST_ROOT' ",
      "for more information.")

### Returns the path to the IgBLAST installation used by igblastr. In case
### of an external installation, the path is returned with the "obtained_via"
### attribute on it. See .get_external_igblast_root() above in this file for
### more information.
### Checks the returned installation only if it's an "external" one.
### Exported!
get_igblast_root <- function()
{
    ## 1. Look for an "internal" IgBLAST installation that is in use, and
    ##    return it if any.
    igblast_root <- get_internal_igblast_root()
    if (!is.na(igblast_root))
        return(igblast_root)  # naked path

    ## 2. Look for an "external" IgBLAST installation and return it if any.
    igblast_root <- .get_external_igblast_root()
    if (!is.na(igblast_root))
        return(igblast_root)  # path has "obtained_via" attribute

    ## 3. Look for "internal" IgBLAST installations. If there are any, put
    ##    the first one in use (that's the one with the highest version)
    ##    and return it.
    all_versions <- .list_installed_igblast_versions()
    if (length(all_versions) != 0L) {
        version <- all_versions[[1L]]  # highest version
        return(set_internal_igblast_root(version))
    }

    ## All the above have failed!
    stop("No IgBLAST installation found.\n  ",
         wmsg(.what_to_do_if_no_internal_installation_yet()))
}

### Sets the path to the IgBLAST installation to use and returns it.
### This can be either an "internal" or an "external" installation.
### In the former case, 'version_or_path' should be the version of an
### existing internal installation.
### In the latter case, it should be the full path (absolute or relative)
### to the root directory of a valid external installation.
### Always checks the returned installation.
### Exported!
set_igblast_root <- function(version_or_path)
{
    if (!isSingleNonWhiteString(version_or_path))
        stop(wmsg("'version_or_path' must be a single (non-empty) string"))

    version <- numeric_version(version_or_path, strict=FALSE)
    if (!is.na(version)) {
        ## 'version_or_path' is a syntactically valid version number expected
        ## to be that of an existing "internal" IgBLAST installation.
        all_versions <- .list_installed_igblast_versions()
        if (length(all_versions) == 0L)
            stop(wmsg("You don't have any igblastr-managed IgBLAST ",
                      "installation yet. ",
                      .what_to_do_if_no_internal_installation_yet()))
        if (!(version_or_path %in% all_versions)) {
            err_msg1 <- c("The supplied version ('", version_or_path, "') ",
                          "is not the version of an igblastr-managed ",
                          "installation of IgBLAST.")
            all_in_1string <- paste0("\"", all_versions, "\"", collapse=", ")
            err_msg2 <- c("List of igblastr-managed installations of IgBLAST ",
                          "(from newest to oldest version): ",
                          all_in_1string, ".")
            stop(wmsg(err_msg1), "\n  ", wmsg(err_msg2))
        }
        igblast_root <- set_internal_igblast_root(version_or_path)
        return(invisible(igblast_root))
    }
    if (dir.exists(version_or_path)) {
        ## 'version_or_path' is the path to an existing directory
        ## expected to contain a valid "external" IgBLAST installation.
        igblast_root <- .set_external_igblast_root(version_or_path)
        ## Remove the USING file if any.
        internal_roots <- get_internal_igblast_roots()
        if (dir.exists(internal_roots)) {
            using_path <- file.path(internal_roots, "USING")
            if (file.exists(using_path))
                unlink(using_path)
        }
        return(invisible(igblast_root))
    }
    stop(wmsg("'version_or_path' must be either a version number ",
              "(e.g. \"1.22.0\") or the path to an existing directory"))
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


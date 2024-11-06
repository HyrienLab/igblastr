### =========================================================================
### install_igblast()
### -------------------------------------------------------------------------


.IGBLAST_ALL_RELEASES_FTP_DIR <-
    "ftp.ncbi.nih.gov/blast/executables/igblast/release/"


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A small set of low-level utils to help find the precompiled IgBLAST
### on NCBI FTP site, for a given IgBLAST release and OS/arch
###

.get_all_releases <- function()
{
    all_releases_ftp_dir <- .IGBLAST_ALL_RELEASES_FTP_DIR
    listing <- try(suppressWarnings(list_ftp_dir(all_releases_ftp_dir)),
                   silent=TRUE)
    if (inherits(listing, "try-error"))
        stop(wmsg("Cannot open URL '", all_releases_ftp_dir, "'. ",
                  "Are you connected to the internet?"))
    ans <- grep("^[0-9]", listing, value=TRUE)
    ans <- as.character(sort(numeric_version(ans), decreasing=TRUE))
    if ("LATEST" %in% listing)
        ans <- c("LATEST", ans)
    ans
}

.get_release_ftp_dir <- function(release="LATEST")
{
    all_releases <- .get_all_releases()
    if (!isSingleNonWhiteString(release))
        stop(wmsg("'release' must be a single (non-empty) string"))
    if (!(release %in% all_releases)) {
        all_in_1string <- paste0("\"", all_releases, "\"", collapse=", ")
        stop(wmsg("'release' must be \"LATEST\" (recommended), or ",
                  "one of the IgBLAST release versions listed at ",
                  .IGBLAST_ALL_RELEASES_FTP_DIR, " (e.g. \"1.21.0\")."),
             "\n  ",
             wmsg("All available releases: ", all_in_1string, "."),
             "\n  ",
             wmsg("Note that old versions have not been tested and ",
                  "are not guaranteed to be compatible with the ",
                  "igblastr package."))
    }
    paste0(.IGBLAST_ALL_RELEASES_FTP_DIR, release, "/")
}

### Returns the suffix of the precompiled NCBI IgBLAST for the specified
### OS/arch, if that OS/arch is supported. Otherwise, returns NA_character_.
.infer_precompiled_suffix_from_OS_arch <- function(OS_arch)
{
    OS <- tolower(OS_arch[["OS"]])
    arch <- tolower(OS_arch[["arch"]])
    if (grepl("^x86.64$", arch)) {
        if (OS == "linux")   return("-x64-linux.tar.gz")
        if (OS == "windows") return("-x64-win64.tar.gz")
        if (OS == "darwin")  return("-x64-macosx.tar.gz")
    } else if (arch == "arm64") {
        if (OS == "darwin")  return(".dmg")
    }
    NA_character_
}

.stop_on_no_precompiled_ncbi_igblast <- function(ftp_dir, OS_arch)
{
    fmt <- paste0("No pre-compiled IgBLAST program available ",
                  "at ", ftp_dir, " for %s.")
    if (anyNA(OS_arch)) {
        err_msg1 <- sprintf(fmt, "your OS/arch")
    } else {
        err_msg1 <- sprintf(fmt, paste(OS_arch, collapse="-"))
    }
    err_msg2 <- c("If there's an existing IgBLAST installation on your ",
                  "machine, please set environment variable IGBLAST_ROOT ",
                  "to point to it. Otherwise, please perform your own ",
                  "installation of IgBLAST e.g. by downloading and compiling ",
                  "the latest source tarball from ", ftp_dir, " (note that ",
                  "compilation can take between 45 min and 1 hour!), then ",
                  "point IGBLAST_ROOT to it. See '?IGBLAST_ROOT' for more ",
                  "information.")
    stop(wmsg(err_msg1), "\n  ", wmsg(err_msg2))
}

.get_precompiled_ncbi_igblast_name <- function(ftp_dir, OS_arch)
{
    suffix <- .infer_precompiled_suffix_from_OS_arch(OS_arch)
    if (is.na(suffix))
        .stop_on_no_precompiled_ncbi_igblast(ftp_dir, OS_arch)
    pattern <- paste0(PRECOMPILED_NCBI_IGBLAST_PREFIX, ".*", suffix)
    listing <- try(suppressWarnings(list_ftp_dir(ftp_dir)), silent=TRUE)
    if (inherits(listing, "try-error"))
        stop(wmsg("Cannot open URL '", ftp_dir, "'. ",
                  "Are you connected to the internet?"))
    idx <- grep(paste0("^", pattern, "$"), listing)
    if (length(idx) == 0L)
        stop(wmsg("Anomaly: no file with a name matching ",
                  pattern, " found at ", ftp_dir))
    listing[[idx[[1L]]]]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_igblast()
###

.projected_igblast_root <- function(ncbi_igblast_name)
{
    internal_roots <- get_internal_igblast_roots()
    version <- infer_igblast_version_from_ncbi_name(ncbi_igblast_name)
    file.path(internal_roots, version)
}

.stop_on_existing_installation <- function(release, igblast_root)
{
    if (release == "LATEST") {
        what <- "LATEST IgBLAST"
    } else {
        what <- c("IgBLAST ", release)
    }
    msg <- c(what, " is already installed in ", igblast_root, "/")
    what_to_do <- " 'force=TRUE' to reinstall."
    internal_igblast_root <- get_internal_igblast_root()
    is_already_in_use <- !is.na(internal_igblast_root) &&
                         igblast_root == internal_igblast_root
    if (is_already_in_use) {
        what_to_do <- c("Use", what_to_do)
    } else {
        version <- basename(igblast_root)
        what_to_do <- c("Call 'set_igblast_root(\"", version, "\")' ",
                        "to use this installation (see '?set_igblast_root' ",
                        "for the details), or use", what_to_do)
    }
    stop(wmsg(msg), "\n  ", wmsg(what_to_do))
}

### Returns the version of IgBLAST being installed which is also the
### basename of its installation directory.
.extract_to_internal_roots <- function(downloaded_file, ncbi_igblast_name)
{
    ## Create "internal roots" folder if it doesn't exist yet.
    internal_roots <- get_internal_igblast_roots()
    if (!dir.exists(internal_roots))
        dir.create(internal_roots, recursive=TRUE)
    if (has_suffix(ncbi_igblast_name, ".tar.gz")) {
        extract_igblast_tarball(downloaded_file, ncbi_igblast_name,
                            destdir=internal_roots)
    } else if (has_suffix(ncbi_igblast_name, ".dmg")) {
        extract_igblast_dmg(downloaded_file, ncbi_igblast_name,
                            destdir=internal_roots)
    } else {
        stop(wmsg("Anomaly: file to extract must be ",
                  "a tarball (.tar.gz file) or .dmg file"))
    }
    infer_igblast_version_from_ncbi_name(ncbi_igblast_name)
}

### TODO: Bad things will happen if more than one R process run
### install_igblast() concurrently. Standard way to address this
### is to put a lock on the "internal roots" folder to get
### exclusive write access to it for the duration of the
### .extract_to_internal_roots() and set_internal_igblast_root() steps.
install_igblast <- function(release="LATEST", force=FALSE, ...)
{
    ftp_dir <- .get_release_ftp_dir(release)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    OS_arch <- get_OS_arch()
    ncbi_igblast_name <- .get_precompiled_ncbi_igblast_name(ftp_dir, OS_arch)
    proj_igblast_root <- .projected_igblast_root(ncbi_igblast_name)
    if (dir.exists(proj_igblast_root)) {
        if (!force)
            .stop_on_existing_installation(release, proj_igblast_root)
    }

    downloaded_file <- download_ftp_file(ftp_dir, ncbi_igblast_name, ...)

    ## Note that bad things will happen if another R process is running
    ## this step at the same time!
    message("Installing IgBLAST at ", proj_igblast_root, " ... ",
            appendLF=FALSE)
    version <- .extract_to_internal_roots(downloaded_file, ncbi_igblast_name)
    message("ok")

    ## Note that bad things will happen if another R process is running
    ## this step at the same time!
    igblast_root <- set_internal_igblast_root(version)

    stopifnot(identical(igblast_root, proj_igblast_root))  # sanity check
    message("IgBLAST ", version, " successfully installed.")
    invisible(igblast_root)
}


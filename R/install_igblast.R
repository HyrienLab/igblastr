### =========================================================================
### install_igblast()
### -------------------------------------------------------------------------


.IGBLAST_ALL_RELEASES_FTP_DIR <-
    "ftp.ncbi.nih.gov/blast/executables/igblast/release/"

.IGBLAST_TARBALL_NAME_PREFIX <- "ncbi-igblast-"

.get_igblast_tarball_name_pattern <- function()
    sprintf("^(%s([.0-9]+)).*$", .IGBLAST_TARBALL_NAME_PREFIX)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A small set of low-level utils to help find the precompiled IgBlast
### tarball on NCBI FTP site, for a given IgBlast release and platform
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
                  "one of the IgBlast release versions listed at ",
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

.how_to_install_manually <- function(what)
{
    msg <- c("Please download and install ", what, " manually from ",
             .get_release_ftp_dir(), ", and set environment ",
             "variable IGBLAST_ROOT accordingly. See '?IGBLAST_ROOT' ",
             "for more information.")
    paste(msg, collapse="")
}

.infer_precompiled_tarball_suffix_from_OS_arch <- function(OS_arch, ftp_dir)
{
    fmt <- paste0("No pre-compiled IgBlast tarball (.tar.gz) ",
                  "available at ", ftp_dir, " for %s.")
    if (anyNA(OS_arch))
        stop(wmsg(sprintf(fmt, "your OS/arch")))
    OS_arch <- paste(OS_arch, collapse="-")
    err_msg <- sprintf(fmt, OS_arch)
    switch(OS_arch,
        `Linux-x86_64`="x64-linux.tar.gz",
        `Windows-x86-64`="x64-win64.tar.gz",
        `Darwin-x86_64`="x64-macosx.tar.gz",
        `Darwin-arm64`=stop(wmsg(
            err_msg, " ",
            .how_to_install_manually("ncbi-igblast-X.Y.Z+.dmg")
        )),
        stop(wmsg(err_msg))
    )
}

.get_precompiled_tarball_name <- function(ftp_dir, OS_arch)
{
    suffix <- .infer_precompiled_tarball_suffix_from_OS_arch(OS_arch, ftp_dir)
    pattern <- paste0(.IGBLAST_TARBALL_NAME_PREFIX, ".*-", suffix)
    listing <- try(suppressWarnings(list_ftp_dir(ftp_dir)), silent=TRUE)
    if (inherits(listing, "try-error"))
        stop(wmsg("Cannot open URL '", ftp_dir, "'. ",
                  "Are you connected to the internet?"))
    idx <- grep(paste0("^", pattern, "$"), listing)
    if (length(idx) == 0L)
        stop(wmsg("Anomaly: no tarball with a name matching ",
                  pattern, " found at ", ftp_dir))
    listing[[idx[[1L]]]]
}

.infer_rootbasename_from_igblast_tarball_name <- function(tarball_name)
{
    pattern <- .get_igblast_tarball_name_pattern()
    sub(pattern, "\\1", tarball_name)
}

.infer_version_from_igblast_tarball_name <- function(tarball_name)
{
    pattern <- .get_igblast_tarball_name_pattern()
    sub(pattern, "\\2", tarball_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_igblast()
###

.projected_igblast_root <- function(tarball_name)
{
    local_executables_dir <- igblastr_local_executables_dir()
    rootbasename <- .infer_rootbasename_from_igblast_tarball_name(tarball_name)
    file.path(local_executables_dir, rootbasename)
}

.stop_on_existing_installation <- function(release, igblast_root)
{
    if (release == "LATEST") {
        what <- "LATEST IgBlast"
    } else {
        what <- c("IgBlast ", release)
    }
    msg <- c(what, " is already installed in ", igblast_root, "/")
    what_to_do <- " 'force=TRUE' to reinstall."
    internal_igblast_root <- get_internal_igblast_root()
    is_already_in_use <- !is.na(internal_igblast_root) &&
                         igblast_root == internal_igblast_root
    if (is_already_in_use) {
        what_to_do <- c("Use", what_to_do)
    } else {
        rootbasename <- basename(igblast_root)
        what_to_do <- c("Call 'set_igblast_root(\"", rootbasename, "\")' ",
                        "to use this installation (see '?set_igblast_root' ",
                        "for the details), or use", what_to_do)
    }
    stop(wmsg(msg), "\n  ", wmsg(what_to_do))
}

.download_tarball <- function(ftp_dir, tarball_file, ...)
{
    url <- paste0(ftp_dir, tarball_file)
    destfile <- tempfile(fileext=".tar.gz")
    code <- try(suppressWarnings(download.file(url, destfile, ...)),
                silent=TRUE)
    if (inherits(code, "try-error") || code != 0L)
        stop(wmsg("Failed to download ", tarball_file, " ",
                  "from ", ftp_dir, ". ",
                  "Are you connected to the internet?"))
    destfile
}

.extract_to_local_executables_dir <- function(tarfile, tarball_name)
{
    local_executables_dir <- igblastr_local_executables_dir()
    ## untar() will create 'local_executables_dir' if it doesn't exist yet.
    code <- suppressWarnings(untar(tarfile, exdir=local_executables_dir))
    if (code != 0L)
        stop(wmsg("Anomaly: something went wrong during ",
                  "extraction of '", tarfile, "' (the local copy of ",
                  "'", tarball_name, "') to '", local_executables_dir, "'."))
    .infer_rootbasename_from_igblast_tarball_name(tarball_name)
}

### TODO: Bad things will happen if more than one R process run
### install_igblast() concurrently. Standard way to address this
### is to put a lock on the igblastr_local_executables_dir() folder
### to get exclusive write access to it for the duration of the
### .extract_to_local_executables_dir() and set_internal_igblast_root() steps.
install_igblast <- function(release="LATEST", force=FALSE, ...)
{
    ftp_dir <- .get_release_ftp_dir(release)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    OS_arch <- get_OS_arch()
    tarball_name <- .get_precompiled_tarball_name(ftp_dir, OS_arch)
    proj_igblast_root <- .projected_igblast_root(tarball_name)
    if (dir.exists(proj_igblast_root)) {
        if (!force)
            .stop_on_existing_installation(release, proj_igblast_root)
        ## Do NOT nuke the existing installation if 'force' is TRUE.
        ## untar() should be able to overwrite it anyways. A major
        ## inconvenient of nuking it now is that if something goes
        ## wrong during .download_tarball() then we'll be left with
        ## nothing. Also if we nuke it now we should also remove
        ## <local_executables_dir>/USING otherwise it will be pointing
        ## to a non-existing installation if download fails.
        ## If we really want to nuke <proj_igblast_root>, we should do
        ## it in .extract_to_local_executables_dir() right before calling
        ## untar(). Or, even safer, we should move <proj_igblast_root>
        ## to a temporary place and nuke it if the new installation is
        ## successful or restore it if it's not. Same mechanism as with
        ## package installation in R.
    }
    tmptarfile <- .download_tarball(ftp_dir, tarball_name, ...)
    ## Note that bad things will happen if another R process is running
    ## the two steps below at the same time!
    rootbasename <- .extract_to_local_executables_dir(tmptarfile, tarball_name)
    igblast_root <- set_internal_igblast_root(rootbasename)
    stopifnot(identical(igblast_root, proj_igblast_root))  # sanity check
    message("IgBlast successfully installed at '", igblast_root, "'.")
    invisible(igblast_root)
}


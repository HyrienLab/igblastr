### =========================================================================
### install_igblast()
### -------------------------------------------------------------------------


.IGBLAST_RELEASE_FTP_DIR <-
    "ftp.ncbi.nih.gov/blast/executables/igblast/release/"

### Returns a single string in the "<OS name>-<arch>" format, or an
### NA_character_ if this information is not available.
### For example, it will return:
###   - "Linux-x86_64" on Intel Linux,
###   - "windows-x64" on Intel Windows,
###   - "Darwin-x86_64" on Intel Mac,
###   - "Darwin-arm64" on ARM Mac (Mac Silicon).
.get_platform <- function()
{
    sys_info <- Sys.info()
    sysname <- sys_info[["sysname"]]
    machine <- sys_info[["machine"]]
    if (!isSingleNonWhiteString(sysname) || !isSingleNonWhiteString(machine))
        return(NA_character_)
    paste0(sysname, "-", machine)
}

.infer_precompiled_tarball_suffix_from_platform <- function(platform)
{
    fmt <- paste0("No precompiled IgBlast tarball available at ",
                  .IGBLAST_RELEASE_FTP_DIR, " for %s.")
    if (is.na(platform))
        stop(wmsg(sprintf(fmt, "your OS/arch")))
    err_msg <- sprintf(fmt, platform)
    switch(platform,
        `Linux-x86_64`="x64-linux.tar.gz",
        `windows-x64`="x64-win64.tar.gz",
        `Darwin-x86_64`="x64-macosx.tar.gz",
        `Darwin-arm64`=stop(wmsg(c(err_msg, " Please download and install ",
                                   "ncbi-igblast-1.22.0+.dmg manually and ",
                                   "set environment variable IGBLAST_ROOT ",
                                   "accordingly. See '?IGBLAST_ROOT' for ",
                                   "more information."))),
        stop(wmsg(err_msg))
    )
}

.get_precompiled_tarball_name <- function(ftp_dir, platform)
{
    suffix <- .infer_precompiled_tarball_suffix_from_platform(platform)
    pattern <- paste0("ncbi-igblast-.*-", suffix)
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

.download_tarball <- function(ftp_dir, tarball_file, ...)
{
    url <- paste0(ftp_dir, tarball_file)
    destfile <- tempfile(fileext=".tar.gz")
    code <- try(suppressWarnings(download.file(url, destfile, ...)),
                silent=TRUE)
    if (inherits(code, "try-error") || code != 0L)
        stop(wmsg("Failed to download '", tarball_file, "' ",
                  "from '", ftp_dir, "'. ",
                  "Are you connected to the internet?"))
    destfile
}

.extract_local_tarball <- function(local_tarball, tarball_name)
{
    exdir <- R_user_dir("igblastr", "cache")
    code <- suppressWarnings(untar(local_tarball, exdir=exdir))
    if (code != 0L)
        stop(wmsg("Anomaly: something went wrong during ",
                  "extraction of '", local_tarball, "', the local copy ",
                  "of '", tarball_name, "', to '", exdir, "'."))
    exdir
}

.infer_igblast_root_basename_from_tarball_name <- function(tarball_name)
{
    pattern <- "^(ncbi-igblast-[^-]*)-.*$"
    sub(pattern, "\\1", tarball_name)
}

.check_igblast_root <- function(igblast_root)
{
    if (!dir.exists(igblast_root))
        stop(wmsg("Anomaly: no IgBlast installation ",
                  "at '", igblast_root, "' (directory doesn't exist)"))
    bin_dir <- file.path(igblast_root, "bin")
    if (!dir.exists(bin_dir))
        stop(wmsg("Anomaly: invalid IgBlast installation ",
                  "at '", igblast_root, "' ",
                  "(no 'bin' subdirectory)"))
    bin_files <- list.files(bin_dir)
    required_bin_files <- c("igblastn", "igblastp",
                            "makeblastdb", "edit_imgt_file.pl")
    for (file in required_bin_files) {
        if (!(file %in% bin_files))
            stop(wmsg("Anomaly: invalid IgBlast installation ",
                      "at '", igblast_root, "' ",
                      "(no '", file, "' file in 'bin' subdirectory)"))
    }
    ## Check that the 'igblastn' and 'igblastp' executables work.
    igblastn_exe <- file.path(bin_dir, "igblastn")
    if (!system_command_works(igblastn_exe, "-version"))
        stop(wmsg("Anomaly: invalid IgBlast installation ",
                  "at '", igblast_root, "' ",
                  "('", igblastn_exe, " -version' does not work)"))
    igblastp_exe <- file.path(bin_dir, "igblastp")
    if (!system_command_works(igblastp_exe, "-version"))
        stop(wmsg("Anomaly: invalid IgBlast installation ",
                  "at '", igblast_root, "' ",
                  "('", igblastp_exe, " -version' does not work)"))
    invisible(igblast_root)
}

install_igblast <- function(release="LATEST", ...)
{
    if (!isSingleNonWhiteString(release))
        stop(wmsg("'release' must be a single string indicating ",
                  "which IgBlast release to install, ",
                  "e.g. \"LATEST\" (recommended) or \"1.22.0\". ",
                  "See ", .IGBLAST_RELEASE_FTP_DIR, " for the list ",
                  "of available releases."))
    platform <- .get_platform()
    ftp_dir <- paste0(.IGBLAST_RELEASE_FTP_DIR, release, "/")
    tarball_name <- .get_precompiled_tarball_name(ftp_dir, platform)
    localfile <- .download_tarball(ftp_dir, tarball_name, ...)
    root_dirname <- .extract_local_tarball(localfile, tarball_name)
    root_basename <-
        .infer_igblast_root_basename_from_tarball_name(tarball_name)
    igblast_root <- file.path(root_dirname, root_basename)
    .check_igblast_root(igblast_root)
    message("IgBlast successfully installed at '", igblast_root, "'.")
    invisible(igblast_root)
}


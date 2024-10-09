.IMGT_VQUEST_REFERENCE_DIRECTORY <-
    "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"

.IG_FILES <- paste0("IG",
    c("HV", "HD", "HJ", "KV", "KJ", "LV", "LJ"), ".fasta")

.TR_FILES <- paste0("TR",
    c("AV", "AJ", "BV", "BD", "BJ", "DV", "DD", "DJ", "GV", "GJ"), ".fasta")

.url_exists <- function(url)
{
    response <- try(HEAD(url), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response), "  Please check your internet connection.")
    response$status_code != 404L
}

normalize_IMGT_VQUEST_organism <- function(organism)
{
    if (!isSingleString(organism))
        stop(wmsg("'organism' must be a single string"))
    if (grepl("^\\s*$", organism))
        stop(wmsg("'organism' contains only whitespace characters"))
    chartr(" ", "_", organism)
}

.get_IMGT_VQUEST_orgdir_url <- function(organism)
{
    organism <- normalize_IMGT_VQUEST_organism(organism)
    orgdir_url <- paste0(.IMGT_VQUEST_REFERENCE_DIRECTORY, organism, "/")
    if (!.url_exists(orgdir_url))
        stop(organism, ": no such organism\n  ",
             wmsg("See ", .IMGT_VQUEST_REFERENCE_DIRECTORY, " for ",
                  "the list of supported organisms."))
    orgdir_url
}

### 'orgdir_url', 'dry.run', and 'destdir' are trusted.
### In particular:
### - 'orgdir_url' is trusted to have the trailing slash, which
###   it should if it was obtained with .get_IMGT_VQUEST_orgdir_url().
### - 'destdir' is trusted to exist.
.download_IG_or_TR_files <- function(orgdir_url, subdir=c("IG", "TR"),
                                     dry.run=FALSE, destdir=".",
                                     method, quiet=FALSE)
{
    subdir <- match.arg(subdir)
    files <- switch(subdir, "IG"=.IG_FILES, "TR"=.TR_FILES)
    subdir_url <- paste0(orgdir_url, subdir, "/")
    downloaded_files <- character(0)
    for (file in files) {
        file_url <- paste0(subdir_url, file)
        ## Not all organisms have all the IG or TR files.
        if (.url_exists(file_url)) {
            if (!dry.run) {
                destfile <- file.path(destdir, file)
                download.file(file_url, destfile, method, quiet)
            }
            downloaded_files <- c(downloaded_files, file)
        }
    }
    downloaded_files
}

### Thin wrapper to .download_IG_or_TR_files() above that creates the IG
### or TR subdir inside 'destdir' and download the files to it.
### We're trying to achieve atomic behavior: in case of error during download,
### the IG or TR subdir is removed. So we get everything or nothing!
### However, note that in case of user interrupt (CTRL+C), we can end up
### with a partial download and corrupted files!
### FIXME: Make the behavior **really** atomic, under any circumstances.
### Like for .download_IG_or_TR_files() above, 'orgdir_url' and 'destdir'
### are trusted.
.download_IG_or_TR_subdir <- function(orgdir_url, subdir=c("IG", "TR"),
                                      destdir=".", method, quiet=FALSE)
{
    subdir <- match.arg(subdir)
    destsubdir <- file.path(destdir, subdir)
    unlink(destsubdir, recursive=TRUE)
    if (!suppressWarnings(dir.create(destsubdir)))
        stop(wmsg("failed to create '", destsubdir, "'"))
    files <- try(.download_IG_or_TR_files(orgdir_url, subdir=subdir,
                                          destdir=destsubdir,
                                          method=method, quiet=quiet))
    if (inherits(files, "try-error")) {
        unlink(destsubdir, recursive=TRUE)
        stop(as.character(files))
    }
    files
}

download_IMGT_germline_sequences <-
    function(organism="Homo_sapiens", subdir=c("IG", "TR", "both"),
             dry.run=FALSE, destdir=".", method, quiet=FALSE)
{
    orgdir_url <- .get_IMGT_VQUEST_orgdir_url(organism)
    subdir <- match.arg(subdir)
    if (!isTRUEorFALSE(dry.run))
        stop(wmsg("'dry.run' must be TRUE or FALSE"))
    if (!isSingleString(destdir))
        stop(wmsg("'destdir' must be a single string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (missing(method))
        method <- getOption("download.file.method", "auto")

    if (subdir == "both") {
        IG_files <- .download_IG_or_TR_files(orgdir_url, subdir="IG",
                                             dry.run=dry.run, destdir=destdir,
                                             method=method, quiet=quiet)
        TR_files <- .download_IG_or_TR_files(orgdir_url, subdir="TR",
                                             dry.run=dry.run, destdir=destdir,
                                             method=method, quiet=quiet)
        files <- c(IG_files, TR_files)
    } else {
        files <- .download_IG_or_TR_files(orgdir_url, subdir=subdir,
                                          dry.run=dry.run, destdir=destdir,
                                          method=method, quiet=quiet)
    }
    if (dry.run) files else invisible(files)
}

IMGT_VQUEST_orgdir_cache <- function(organism)
{
    organism <- normalize_IMGT_VQUEST_organism(organism)
    file.path(R_user_dir("igblastr", "cache"),
              "IMGT_V-QUEST_reference_directory", organism)
}

download_IMGT_germline_sequences_to_cache <-
    function(organism="Homo_sapiens", subdir=c("IG", "TR", "both"),
             method, quiet=FALSE)
{
    organism <- normalize_IMGT_VQUEST_organism(organism)
    orgdir_url <- .get_IMGT_VQUEST_orgdir_url(organism)
    subdir <- match.arg(subdir)
    if (missing(method))
        method <- getOption("download.file.method", "auto")

    orgdir_cache <- IMGT_VQUEST_orgdir_cache(organism)
    if (!dir.exists(orgdir_cache))
        dir.create(orgdir_cache, recursive=TRUE)

    if (subdir == "both") {
        IG_files <- .download_IG_or_TR_subdir(orgdir_url, subdir="IG",
                                              destdir=orgdir_cache,
                                              method=method, quiet=quiet)
        TR_files <- .download_IG_or_TR_subdir(orgdir_url, subdir="TR",
                                              destdir=orgdir_cache,
                                              method=method, quiet=quiet)
        files <- c(IG_files, TR_files)
    } else {
        files <- .download_IG_or_TR_subdir(orgdir_url, subdir=subdir,
                                           destdir=orgdir_cache,
                                           method=method, quiet=quiet)
    }
    files
}


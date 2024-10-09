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

.get_IMGT_VQUEST_orgdir_url <- function(organism)
{
    if (!isSingleString(organism))
        stop(wmsg("'organism' must be a single string"))
    if (grepl("^\\s*$", organism))
        stop(wmsg("'organism' contains only whitespace characters"))
    organism <- chartr(" ", "_", organism)
    orgdir_url <- paste0(.IMGT_VQUEST_REFERENCE_DIRECTORY, organism, "/")
    if (!.url_exists(orgdir_url))
        stop(organism, ": no such organism\n  ",
             wmsg("See ", .IMGT_VQUEST_REFERENCE_DIRECTORY, " for ",
                  "the list of supported organisms."))
    orgdir_url
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

    if (subdir == "both") {
        downloaded_IG_files <- download_IMGT_germline_sequences(organism,
                                   subdir="IG", dry.run=dry.run,
                                   destdir=destdir, method=method, quiet=quiet)
        downloaded_TR_files <- download_IMGT_germline_sequences(organism,
                                   subdir="TR", dry.run=dry.run,
                                   destdir=destdir, method=method, quiet=quiet)
        downloaded_files <- c(downloaded_IG_files, downloaded_TR_files)
    } else {
        files <- if (subdir == "IG") .IG_FILES else .TR_FILES
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
    }
    if (dry.run) downloaded_files else invisible(downloaded_files)
}


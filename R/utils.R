.IMGT_V_QUEST_REFERENCE_DIRECTORY <-
    "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory"

.IG_FILES <- paste0("IG",
    c("HV", "HD", "HJ", "KV", "KJ", "LV", "LJ"), ".fasta")

.TR_FILES <- paste0("TR",
    c("AV", "AJ", "BV", "BD", "BJ", "DV", "DD", "DJ", "GV", "GV"), ".fasta")

url_exists <- function(url)
{
    code <- HEAD(url)$status_code
    code >= 200L && code < 300L
}

download_IMGT_germline_sequences <-
    function(organism="Homo_sapiens", subdir=c("IG", "TR", "both"),
             destdir=".", method, quiet=FALSE)
{
    if (!isSingleString(organism))
        stop(wmsg("'organism' must be a single string"))
    subdir <- match.arg(subdir)
    if (!isSingleString(destdir))
        stop(wmsg("'destdir' must be a single string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (subdir == "both") {
        downloaded_IG_files <- download_IMGT_germline_sequences(organism,
                                   subdir="IG", destdir=destdir, quiet=quiet)
        downloaded_TR_files <- download_IMGT_germline_sequences(organism,
                                   subdir="TR", destdir=destdir, quiet=quiet)
        downloaded_files <- c(downloaded_IG_files, downloaded_TR_files)
    } else {
        files <- if (subdir == "IG") .IG_FILES else .TR_FILES
        subdir_url <- paste(.IMGT_V_QUEST_REFERENCE_DIRECTORY,
                            organism, subdir, sep="/")
        downloaded_files <- character(0)
        for (file in files) {
            file_url <- paste0(subdir_url, "/", file)
            ## Not all organisms have all the IG or TR files.
            if (url_exists(file_url)) {
                destfile <- file.path(destdir, file)
                download.file(file_url, destfile, method, quiet)
                downloaded_files <- c(downloaded_files, file)
            }
        }
    }
    invisible(downloaded_files)
}


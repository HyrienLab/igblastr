.get_igblast_root <- function()
{
    igblast_root <- Sys.getenv("IGBLAST_ROOT")
    if (!nzchar(igblast_root))
        stop(wmsg("environment variable IGBLAST_ROOT must be set"))
    igblast_root
}

.get_igblastn_exe <- function(igblast_root)
{
    igblastn_exe <- file.path(igblast_root, "bin", "igblastn")
    if (!file.exists(igblastn_exe))
        stop(wmsg("invalid IGBLAST_ROOT"))
    igblastn_exe
}

igblastn <- function(query, args=character())
{
    igblast_root <- .get_igblast_root()
    igblastn_exe <- .get_igblastn_exe(igblast_root)

    if (!isSingleNonWhiteString(query))
        stop(wmsg("'query' must be a single string ",
                  "containing the path to the input file"))
    query <- file_path_as_absolute(query)

    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(igblast_root)

    system2(igblastn_exe, c(paste("-query", query), args))
}



igblastn <- function(query, args=character())
{
    igblast_root <- get_igblast_root()
    igblastn_exe <- get_igblast_exe("igblastn")

    if (!isSingleNonWhiteString(query))
        stop(wmsg("'query' must be a single string ",
                  "containing the path to the input file"))
    query <- file_path_as_absolute(query)

    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(igblast_root)

    system2(igblastn_exe, c(paste("-query", query), args))
}


### =========================================================================
### igblast_info() and related
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### has_igblast()
###

has_igblast <- function()
{
    igblast_root <- try(get_igblast_root(), silent=TRUE)
    !inherits(igblast_root, "try-error")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblast_version()
### list_igblast_internal_data()
### igblast_info()
###

.extract_version_from_cmd_output <- function(output)
{
    sub("^igblastn: *", "", output[[1L]])
}

igblast_version <- function()
{
    igblastn_exe <- get_igblast_exe("igblastn")
    out <- system2(igblastn_exe, "-version", stdout=TRUE)
    .extract_version_from_cmd_output(out)
}

list_igblast_internal_data <- function()
{
    igblast_root <- get_igblast_root()
    internal_data <- file.path(igblast_root, "internal_data")
    if (!dir.exists(internal_data))
        return(character(0))
    list.files(internal_data)
}

print.igblast_info <- function(x, ...)
{
    x <- lapply(x, function(x) paste(x, collapse="; "))
    x <- paste0(names(x), ": ", as.character(x))
    cat(x, sep="\n")
}

igblast_info <- function()
{
    igblast_root <- get_igblast_root()
    OS_arch <- get_OS_arch()
    OS <- OS_arch[["OS"]]
    igblastn_exe <- make_igblast_exe_path(igblast_root, cmd="igblastn", OS=OS)
    igblastn_version <- system2(igblastn_exe, "-version", stdout=TRUE)
    version <- .extract_version_from_cmd_output(igblastn_version)
    #igblastp_exe <- make_igblast_exe_path(igblast_root, cmd="igblastp", OS=OS)
    #igblastp_version <- system2(igblastp_exe, "-version", stdout=TRUE)
    build <- igblastn_version[[2L]]  # should be same as igblastp_version[[2L]]
    build <- sub("^ *Package: ", "", build)
    internal_data <- list_igblast_internal_data()
    if (length(internal_data) == 0L) {
        internal_data <- "none!"
    } else {
        internal_data <- paste0(internal_data, collapse=", ")
    }

    ans <- list(
        igblast_root=igblast_root,
        Version=version,
        `OS/arch`=paste(OS_arch, collapse="/"),
        Build=build,
        #igblastn_exe=igblastn_exe,
        #igblastn_version=igblastn_version[[1L]],
        #igblastp_exe=igblastp_exe,
        #igblastp_version=igblastp_version[[1L]],
        internal_data=internal_data
    )
    class(ans) <- "igblast_info"
    ans
}


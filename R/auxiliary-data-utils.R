### =========================================================================
### Manipulation of IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_igblast_auxiliary_data()
###

get_igblast_auxiliary_data <- function(organism)
{
    organism <- normalize_igblast_organism(organism)
    igblast_root <- get_igblast_root()
    dirpath <- file.path(igblast_root, "optional_file")
    auxiliary_data <- file.path(dirpath, paste0(organism, "_gl.aux"))
    if (!file.exists(auxiliary_data))
        stop(wmsg("no auxiliary data found in ", dirpath, " for ", organism))
    auxiliary_data
}


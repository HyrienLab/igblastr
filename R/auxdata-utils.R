### =========================================================================
### Access IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_germline_db_auxdata_path()
### get_auxdata_path()
###

.get_igblast_auxdata_path <- function(which, organism)
{
    organism <- normalize_igblast_organism(organism)
    auxdata_dir <- file.path(path_to_igdata(which), "optional_file")
    auxdata_filename <- paste0(organism, "_gl.aux")
    auxdata_path <- file.path(auxdata_dir, auxdata_filename)
    if (!file.exists(auxdata_path))
        stop(wmsg("no auxiliary data found in ",
                  auxdata_dir, " for ", organism))
    auxdata_path
}

make_germline_db_auxdata_path <- function(db_path)
{
    stopifnot(dir.exists(db_path))
    ## gl.aux: germline auxiliary file ("gl" in gl.aux stands for germline).
    file.path(db_path, "auxiliary_data", "gl.aux")
}

get_auxdata_path <- function(organism, which=c("live", "original"))
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    which <- match.arg(which)

    if (!valid_germline_db_name(organism))
        return(.get_igblast_auxdata_path(which, organism))

    ## Treat 'organism' as a valid germline db name.
    db_path <- get_germline_db_path(organism)
    auxdata_path <- make_germline_db_auxdata_path(db_path)

    if (!dir.exists(dirname(auxdata_path)))
        stop(wmsg("no auxiliary data found in germline db ", organism))
    if (!file.exists(auxdata_path))
        stop(wmsg("auxiliary data file ", basename(auxdata_path), " ",
                  "not found in germline db ", organism))
    auxdata_path
}

get_igblast_auxiliary_data <- function(...)
{
    .Deprecated("get_auxdata_path")
    get_auxdata_path(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_auxdata()
###

### IMPORTANT NOTE: Unlike with the data.frame returned by load_intdata(),
### all the positions in the data.frame returned by load_auxdata() (that is,
### the positions reported in columns 'coding_frame_start' and 'cdr3_end')
### are 0-based!
load_auxdata <- function(organism, which=c("live", "original"))
{
    which <- match.arg(which)
    auxdata_path <- get_auxdata_path(organism, which=which)
    read_auxdata(auxdata_path)
}

load_igblast_auxiliary_data <- function(...)
{
    .Deprecated("load_auxdata")
    load_auxdata(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_germline_db_auxdata()
###

.do_compute_germline_db_auxdata <- function(db_path, ...)
{
    J_alleles <- readDNAStringSet(get_db_fasta_file(db_path, "J"))
    compute_auxdata(J_alleles, ...)
}

compute_germline_db_auxdata <- function(db_name, ...)
{
    check_germline_db_name(db_name)
    db_path <- get_germline_db_path(db_name)
    .do_compute_germline_db_auxdata(db_path, ...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### add_computed_auxdata_if_missing()
###

### Not exported!
### If the auxdata in germline db at 'db_path' is missing then
### add_computed_auxdata_if_missing() computes and adds it.
### Returns FALSE if the germline db already includes auxdata, and TRUE if
### if doesn't and the auxdata could be computed and added sucessfully.
### Raises an error if the auxdata was missing but the computed auxdata
### could not be added because its "coding_frame_start" column contains NAs.
add_computed_auxdata_if_missing <- function(db_path, verbose=FALSE)
{
    stopifnot(isTRUEorFALSE(verbose))
    auxdata_path <- make_germline_db_auxdata_path(db_path)
    if (file.exists(auxdata_path))
        return(FALSE)
    if (verbose)
        message(wmsg("Computing the auxdata"), " ... ", appendLF=FALSE)
    auxdata <- .do_compute_germline_db_auxdata(db_path)
    ok <- add_computed_auxdata_to_db(auxdata, db_path, verbose=verbose)
    if (!ok)
        stop(wmsg("failed to add the computed auxdata to ", basename(db_path)))
    if (verbose)
        message("ok.")
    TRUE
}


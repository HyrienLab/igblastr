### =========================================================================
### igblastn()
### -------------------------------------------------------------------------


### 'args' must be a named character vector where the names are
### valid 'igblastn' parameter names (e.g. "organism") and the values
### are the parameter values (e.g. "rabbit").
.run_igblastn <- function(igblast_root, args,
                          show.in.browser=FALSE, show.command.only=FALSE)
{
    igblastn_exe <- make_igblast_exe_path(igblast_root, "igblastn")
    if (!is.character(args))
        stop(wmsg("'args' must be character vector"))
    args_names <- names(args)
    if (is.null(args_names))
        stop(wmsg("'args' must have names on it"))
    if (!isTRUEorFALSE(show.in.browser))
        stop(wmsg("'show.in.browser' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.command.only))
        stop(wmsg("'show.command.only' must be TRUE or FALSE"))

    if (show.in.browser) {
        stdout <- tempfile(fileext=".html")
    } else {
        stdout <- ""
    }

    args <- paste0("-", args_names, " ", args)
    if (show.command.only) {
        cmd <- c(igblastn_exe, args)
        cmd_in_1string <- paste(cmd, collapse=" ")
        cat(cmd_in_1string, "\n", file=stdout, sep="")
        if (show.in.browser)
            display_local_file_in_browser(stdout)
        return(invisible(cmd))  # returns the command in a character vector
    }

    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))
    status <- system2(igblastn_exe, args=args, stdout=stdout)
    if (status == 0 && show.in.browser)
        display_local_file_in_browser(stdout)
    invisible(status)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_outfmt()
###

.normarg_outfmt <- function(outfmt="AIRR")
{
    if (isSingleString(outfmt) && toupper(outfmt) == "AIRR")
        return(19)
    if (!(isSingleNumber(outfmt) && outfmt %in% c(3, 4, 7, 19)))
        stop(wmsg("'outfmt' must be \"AIRR\" or one of 3, 4, 7, 19 ",
                  "(19 means \"AIRR\")"))
    as.integer(outfmt)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_organism()
###

### Maps 'db_name' to one of the organism names returned by
### list_igblast_internal_data().
.infer_igblast_organism_from_db_name <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name))
    if (grepl("Homo.sapiens", db_name, ignore.case=TRUE))
        return("human")
    if (grepl("Mus.musculus", db_name, ignore.case=TRUE))
        return("mouse")
    if (grepl("Oryctolagus.cuniculus", db_name, ignore.case=TRUE))
        return("rabbit")
    if (grepl("Rattus.norvegicus", db_name, ignore.case=TRUE))
        return("rat")
    if (grepl("Macaca.mulatta", db_name, ignore.case=TRUE))
        return("rhesus_monkey")
    NA_character_
}

.normarg_organism <- function(organism="auto", db_name)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    if (organism == "auto") {
        organism <- .infer_igblast_organism_from_db_name(db_name)
        if (!is.na(organism))
            return(organism)
        stop(wmsg("Don't know how to infer 'organism' from germline ",
                  "db name \"", db_name, "\". Please set the 'organism' ",
                  "argument to the name of the IgBLAST internal data to use. ",
                  "Use list_igblast_internal_data() to list all valid names."))
    }
    internal_data <- list_igblast_internal_data()
    organism <- tolower(organism)
    if (!(organism %in% list_igblast_internal_data())) {
        all_in_1string <- paste0("\"", internal_data, "\"", collapse=", ")
        stop(wmsg("'organism' must be one of ", all_in_1string))
    }
    organism
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_igblastn_germline_db_args()
###

.make_igblastn_germline_db_args <- function(db_name)
{
    db_path <- get_germline_db_path(db_name)
    VDJ <- c("V", "D", "J")
    setNames(file.path(db_path, VDJ), paste0("germline_db_", VDJ))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn()
###

igblastn <- function(query, outfmt="AIRR", organism="auto", ...,
                     show.in.browser=FALSE, show.command.only=FALSE)
{
    igblast_root <- get_igblast_root()
    db_name <- use_germline_db()
    if (!isSingleNonWhiteString(query))
        stop(wmsg("'query' must be a single string ",
                  "containing the path to the input file"))
    query <- file_path_as_absolute(query)
    outfmt <- .normarg_outfmt(outfmt)
    organism <- .normarg_organism(organism, db_name)
    auxiliary_data <- file.path(igblast_root, "optional_file",
                                paste0(organism, "_gl.aux"))

    ## Turn extra args into a named character vector.
    extra_args <- list(...)
    extra_args <- setNames(as.character(extra_args), names(extra_args))
    germline_db_args <- .make_igblastn_germline_db_args(db_name)

    args <- c(query=query, outfmt=outfmt, organism=organism,
              auxiliary_data=auxiliary_data, extra_args, germline_db_args)

    .run_igblastn(igblast_root, args, show.in.browser=show.in.browser,
                                      show.command.only=show.command.only)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn_help()
###

igblastn_help <- function(long.help=FALSE, show.in.browser=FALSE)
{
    if (!isTRUEorFALSE(long.help))
        stop(wmsg("'long.help' must be TRUE or FALSE"))
    igblast_root <- get_igblast_root()
    args <- if (long.help) c(help="") else c(h="")
    .run_igblastn(igblast_root, args, show.in.browser=show.in.browser)
}


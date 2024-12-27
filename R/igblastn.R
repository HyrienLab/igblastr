### =========================================================================
### igblastn()
### -------------------------------------------------------------------------


### 'args' must be a named character vector where the names are
### valid 'igblastn' parameter names (e.g. "organism") and the values
### are the parameter values (e.g. "rabbit").
.as_igblastn_args <- function(args)
{
    if (!is.character(args))
        stop(wmsg("'args' must be character vector"))
    args_names <- names(args)
    if (is.null(args_names))
        stop(wmsg("'args' must have names on it"))
    paste0("-", args_names, " ", args)
}

.show_igblastn_command <- function(igblastn_exe, args, show.in.browser=FALSE)
{
    cmd <- c(igblastn_exe, args)
    cmd_in_1string <- paste(cmd, collapse=" ")
    outfile <- if (show.in.browser) tempfile() else ""
    cat(cmd_in_1string, "\n", file=outfile, sep="")
    if (show.in.browser)
        display_local_file_in_browser(outfile)
    cmd  # returns the command in a character vector
}

### 'args' must be a named character vector. See .as_igblastn_args() above
### for details.
.run_igblastn <- function(igblast_root, outfmt, args,
                          show.in.browser=FALSE, show.command.only=FALSE)
{
    igblastn_exe <- make_igblast_exe_path(igblast_root, "igblastn")
    args <- c(paste("-outfmt", outfmt), .as_igblastn_args(args))

    if (!isTRUEorFALSE(show.in.browser))
        stop(wmsg("'show.in.browser' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.command.only))
        stop(wmsg("'show.command.only' must be TRUE or FALSE"))

    if (show.command.only) {
        cmd <- .show_igblastn_command(igblastn_exe, args,
                                      show.in.browser=show.in.browser)
        return(invisible(cmd))  # returns the command in a character vector
    }

    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))

    outfile <- tempfile()
    status <- system2(igblastn_exe, args=args, stdout=outfile)
    if (status != 0)
        stop(wmsg("'igblastn' returned an error"))

    if (outfmt != 19) {
        if (show.in.browser)
            return(display_local_file_in_browser(outfile))
        return(cat(readLines(outfile), sep="\n"))
    }

    ## AIRR output format is tabulated.
    df <- read.table(outfile, header=TRUE, sep="\t")
    if (show.in.browser) {
        temp_html <- tempfile(fileext=".html")
        print(xtable(df), type="html", file=temp_html)
        temp_url <- paste0("file://", temp_html)
        browseURL(temp_url)
    }
    tibble(df)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_query()
###

.normarg_query <- function(query)
{
    if (isSingleNonWhiteString(query))
        return(file_path_as_absolute(query))
    if (is(query, "DNAStringSet")) {
        if (is.null(names(query)))
            stop(wmsg("DNAStringSet object 'query' must have names"))
        filepath <- tempfile(fileext=".fasta")
        writeXStringSet(query, filepath)
        return(filepath)
    }
    stop(wmsg("'query' must be a single (non-empty) string ",
              "that contains the path to the input file. ",
              "Alternatively it can be a named DNAStringSet object."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_outfmt()
###

.normarg_outfmt <- function(outfmt="AIRR")
{
    if (isSingleString(outfmt) && toupper(outfmt) == "AIRR")
        return(19L)
    if (!(isSingleNumber(outfmt) && outfmt %in% c(3, 4, 7, 19)))
        stop(wmsg("'outfmt' must be \"AIRR\" or one of 3, 4, 7, 19 ",
                  "(19 means \"AIRR\")"))
    as.integer(outfmt)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_organism()
###

### Maps 'db_name' to one of the organism names returned by
### list_igblast_organisms().
.infer_igblast_organism_from_germline_db_name <- function(db_name)
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

.normarg_organism <- function(organism="auto", germline_db_name)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    if (organism == "auto") {
        organism <-
            .infer_igblast_organism_from_germline_db_name(germline_db_name)
        if (!is.na(organism))
            return(organism)
        stop(wmsg("Don't know how to infer 'organism' from germline ",
                  "db name \"", germline_db_name, "\". Please set ",
                  "the 'organism' argument to the name of the IgBLAST ",
                  "internal data to use. ",
                  "Use list_igblast_organisms() to list all valid names."))
    }
    all_organisms <- list_igblast_organisms()
    organism <- tolower(organism)
    if (!(organism %in% list_igblast_organisms())) {
        all_in_1string <- paste0("\"", all_organisms, "\"", collapse=", ")
        stop(wmsg("'organism' must be one of ", all_in_1string))
    }
    organism
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_igblastn_germline_db_args()
###

.make_igblastn_germline_db_args <- function(germline_db_name)
{
    db_path <- get_germline_db_path(germline_db_name)
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
    germline_db_name <- use_germline_db()  # cannot be ""
    c_region_db_name <- use_c_region_db()  # can be ""
    query <- .normarg_query(query)
    outfmt <- .normarg_outfmt(outfmt)
    organism <- .normarg_organism(organism, germline_db_name)
    auxiliary_data <- file.path(igblast_root, "optional_file",
                                paste0(organism, "_gl.aux"))

    ## Turn extra args into a named character vector.
    extra_args <- list(...)
    extra_args <- setNames(as.character(extra_args), names(extra_args))
    germline_db_args <- .make_igblastn_germline_db_args(germline_db_name)

    args <- c(query=query, organism=organism,
              auxiliary_data=auxiliary_data, extra_args, germline_db_args)
    if (c_region_db_name != "") {
        c_region_db_path <- get_c_region_db_path(c_region_db_name)
        args <- c(args, c_region_db=file.path(c_region_db_path, "C"))
    }

    .run_igblastn(igblast_root, outfmt, args,
                  show.in.browser=show.in.browser,
                  show.command.only=show.command.only)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn_help()
###

igblastn_help <- function(long.help=FALSE, show.in.browser=FALSE)
{
    if (!isTRUEorFALSE(long.help))
        stop(wmsg("'long.help' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.in.browser))
        stop(wmsg("'show.in.browser' must be TRUE or FALSE"))

    igblast_root <- get_igblast_root()
    igblastn_exe <- make_igblast_exe_path(igblast_root, "igblastn")
    args <- if (long.help) "-help" else "-h"

    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))
    outfile <- if (show.in.browser) tempfile() else ""
    status <- system2(igblastn_exe, args=args, stdout=outfile)
    if (status != 0)
        stop(wmsg("'igblastn' returned an error"))
    if (show.in.browser)
        display_local_file_in_browser(outfile)
}


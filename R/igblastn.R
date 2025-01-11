### =========================================================================
### igblastn()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_query()
###

.check_query_seq_lengths <- function(seqlens)
{
    stopifnot(is.integer(seqlens))
    seqids <- names(seqlens)
    stopifnot(!is.null(seqids))
    empty_idx <- which(seqlens == 0L)
    if (length(empty_idx) != 0L) {
        in1string <- paste(seqids[empty_idx], collapse=", ")
        stop(wmsg("the following sequences in 'query' are empty ",
                  "(showing seq ids): ", in1string))
    }
}

.normarg_query <- function(query)
{
    if (isSingleNonWhiteString(query)) {
        path <- file_path_as_absolute(query)
        .check_query_seq_lengths(fasta.seqlengths(path))
    } else if (is(query, "DNAStringSet")) {
        if (is.null(names(query)))
            stop(wmsg("DNAStringSet object 'query' must have names"))
        .check_query_seq_lengths(setNames(width(query), names(query)))
        path <- tempfile(fileext=".fasta")
        writeXStringSet(query, path)
    } else {
        stop(wmsg("'query' must be a single (non-empty) string ",
                  "that contains the path to the input file (FASTA) ",
                  "or a named DNAStringSet object"))
    }
    path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_outfmt()
###

.stop_on_invalid_outfmt <- function()
{
    msg1 <- c("'outfmt' must be one of 3, 4, 7, 19, \"AIRR\" (\"AIRR\" ",
              "is an alias for 19), or a string describing a customized ",
              "format 7.")
    msg2 <- c("The string describing a customized format 7 must start ",
              "with a 7 followed by the desired hit table fields ",
              "(a.k.a. format specifiers) separated with spaces. ",
              "For example: \"7 std qseq sseq btop\". Note that 'std' ",
              "stands for 'qseqid sseqid pident length mismatch gapopen ",
              "gaps qstart qend sstart send evalue bitscore', which is ",
              "the default.")
    msg3 <- c("Use list_supported_format_specifiers() to list all supported ",
              "format specifiers.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### 'outfmt' is assumed to be a single string.
.check_customized_format_7 <- function(outfmt)
{
    if (!has_prefix(outfmt, "7 "))
        .stop_on_invalid_outfmt()
    user_specifiers <- substr(outfmt, 3L, nchar(outfmt))
    user_specifiers <- strsplit(user_specifiers, " ", fixed=TRUE)[[1L]]
    user_specifiers <- user_specifiers[nzchar(user_specifiers)]
    supported_specifiers <- c("std", names(list_supported_format_specifiers()))
    invalid_specifiers <- setdiff(user_specifiers, supported_specifiers)
    if (length(invalid_specifiers) != 0L) {
        in1string <- paste(invalid_specifiers, collapse=", ")
        stop(wmsg("invalid format specifier(s): ", in1string))
    }
}

### 'outfmt' can be 3, 4, 7, 19, "AIRR", or a single string describing a
### customized format 7 e.g. "7 std qseq sseq btop".
### See .stop_on_invalid_outfmt() above for more information.
### Returns a single string.
.normarg_outfmt <- function(outfmt=7)
{
    if (isSingleNumber(outfmt)) {
        if (!(outfmt %in% c(3, 4, 7, 19)))
            .stop_on_invalid_outfmt()
        outfmt <- as.character(as.integer(outfmt))
    } else if (isSingleNonWhiteString(outfmt)) {
        outfmt <- trimws(outfmt)
        if (toupper(outfmt) == "AIRR") {
            outfmt <- "19"
        } else {
            if (!(outfmt %in% c("3", "4", "7", "19")))
                .check_customized_format_7(outfmt)
        }
    } else {
        .stop_on_invalid_outfmt()
    }
    outfmt
}

### Returns 3, 4, 7, or 19.
.extract_fmt_nb <- function(outfmt)
{
    stopifnot(isSingleNonWhiteString(outfmt))
    fmt_nb <- strsplit(trimws(outfmt), " ", fixed=TRUE)[[1L]][[1L]]
    stopifnot(fmt_nb %in% c("3", "4", "7", "19"))
    as.integer(fmt_nb)
}

print.igblastn_raw_output <- function(x, ...) cat(x, sep="\n")


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
### .normarg_seqidlist()
###

.normarg_seqidlist_file <- function(seqidlist)
{
    stopifnot(inherits(seqidlist, "file"))
    path <- summary(seqidlist)$description
    ## Check that the file exists and is readable.
    res <- try(suppressWarnings(open(seqidlist, open="r")), silent=TRUE)
    if (inherits(res, "try-error"))
        stop(wmsg("unable to open file '", path, "' for reading"))
    close(seqidlist)
    path
}

.normarg_seqidlist_character <- function(seqidlist)
{
    stopifnot(is.character(seqidlist))
    if (anyNA(seqidlist))
        stop(wmsg("character vector 'seqidlist' contains NAs"))
    seqidlist <- trimws(seqidlist)
    if (!all(nzchar(seqidlist)))
        stop(wmsg("character vector 'seqidlist' contains ",
                  "empty or white strings"))
    path <- tempfile()
    writeLines(seqidlist, path)
    path
}

.check_user_seqids <- function(user_seqids, region_type)
{
    sequences <- load_germline_db(use_germline_db(), region_types=region_type)
    valid_seqids <- names(sequences)
    invalid_seqids <- setdiff(user_seqids, valid_seqids)
    if (length(invalid_seqids) != 0L) {
        in1string <- paste0(invalid_seqids, collapse=", ")
        msg1 <- c("Seq id(s) not found in germline ",
                  "database ", region_type, ": ", in1string)
        msg2a <- "Note that you can use:"
        code <- c("names(load_germline_db(use_germline_db(), ",
                  "region_types=\"", region_type, "\"))")
        msg2b <- c("to obtain the list of valid seq ids ",
                   "from germline database ", region_type, ".")
        stop(wmsg(msg1), "\n  ",
             wmsg(msg2a), "\n    ", code, "\n  ", wmsg(msg2b))
    }
}

.normarg_seqidlist <- function(seqidlist, region_type=c("V", "D", "J"))
{
    region_type <- match.arg(region_type)
    if (is.null(seqidlist))
        return(NULL)
    if (inherits(seqidlist, "file")) {
        path <- .normarg_seqidlist_file(seqidlist)
    } else if (is.character(seqidlist)) {
        path <- .normarg_seqidlist_character(seqidlist)
    } else {
        stop(wmsg("'seqidlist' must be NULL, or a 'file' object, ",
                  "or a character vector"))
    }
    user_seqids <- readLines(path)
    .check_user_seqids(user_seqids, region_type)
    path
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
### .parse_igblastn_output()
###

### TODO: Parse output format 3 and 4.
.parse_igblastn_output <- function(out_file, fmt_nb)
{
    out <- readLines(out_file)
    if (fmt_nb == 7)
        return(parse_fmt7(out))
    warning("parsing of igblastn output format ", fmt_nb, " is not ready yet")
    class(out) <- "igblastn_raw_output"
    out
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn()
###

### 'args' must be a named character vector where the names are
### valid 'igblastn' parameter names (e.g. "organism") and the values
### are the parameter values (e.g. "rabbit").
.make_command_line_args <- function(args)
{
    if (!is.character(args))
        stop(wmsg("'args' must be character vector"))
    args_names <- names(args)
    if (is.null(args_names))
        stop(wmsg("'args' must have names on it"))
    quoteme_idx <- grep(" ", args, fixed=TRUE)
    if (length(quoteme_idx) != 0L)
        args[quoteme_idx] <- paste0("'", args[quoteme_idx], "'")
    paste0("-", args_names, " ", args)
}

.show_igblastn_command <- function(igblast_root, cmd_args,
                                   show.in.browser=FALSE)
{
    igblastn_exe <- make_igblast_exe_path(igblast_root, "igblastn")
    cmd <- c(igblastn_exe, cmd_args)
    cmd_in_1string <- paste(cmd, collapse=" ")
    outfile <- if (show.in.browser) tempfile() else ""
    cat(cmd_in_1string, "\n", file=outfile, sep="")
    if (show.in.browser)
        display_local_file_in_browser(outfile)
    cmd  # returns the command in a character vector
}

.parse_warnings_or_errors <- function(lines, prefix)
{
    stopifnot(is.character(lines), isSingleNonWhiteString(prefix))
    keep_idx <- which(has_prefix(tolower(lines), tolower(prefix)))
    lines <- lines[keep_idx]
    msgs <- trimws(substr(lines, nchar(prefix)+1L, nchar(lines)))
    msgs[nzchar(msgs)]
}

.parse_and_issue_warnings <- function(stderr_file)
{
    warn_msgs <- .parse_warnings_or_errors(readLines(stderr_file), "Warning:")
    for (msg in warn_msgs)
        warning(wmsg(msg))
}

.stop_on_igblastn_exe_error <- function(stderr_file)
{
    err_msgs <- .parse_warnings_or_errors(readLines(stderr_file), "Error:")
    if (length(err_msgs) == 0L)  # could this ever happen?
        stop(wmsg("'igblastn' returned an unknown error"))
    err_msgs <- vapply(err_msgs, wmsg, character(1), USE.NAMES=FALSE)
    stop(paste(err_msgs, collapse="\n  "))
}

### 'cmd_args' must be a named character vector. See .make_command_line_args()
### above for details.
### Returns the path to the output file.
.run_igblastn <- function(igblast_root, cmd_args)
{
    igblastn_exe <- make_igblast_exe_path(igblast_root, "igblastn")
    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))

    out_file <- tempfile()
    stderr_file <- tempfile()
    cmd_args <- c(cmd_args, paste("-out", out_file))
    status <- system2(igblastn_exe, args=cmd_args, stderr=stderr_file)
    .parse_and_issue_warnings(stderr_file)
    if (status != 0)
        .stop_on_igblastn_exe_error(stderr_file)
    out_file
}

igblastn <- function(query, outfmt=7, parse.out=TRUE,
                     organism="auto",
                     germline_db_V_seqidlist=NULL,
                     germline_db_D_seqidlist=NULL,
                     germline_db_J_seqidlist=NULL,
                     ...,
                     show.in.browser=FALSE, show.command.only=FALSE)
{
    if (!isTRUEorFALSE(parse.out))
        stop(wmsg("'parse.out' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.in.browser))
        stop(wmsg("'show.in.browser' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.command.only))
        stop(wmsg("'show.command.only' must be TRUE or FALSE"))

    igblast_root <- get_igblast_root()
    germline_db_name <- use_germline_db()  # cannot be ""
    c_region_db_name <- use_c_region_db()  # can be ""
    query <- .normarg_query(query)
    outfmt <- .normarg_outfmt(outfmt)
    fmt_nb <- .extract_fmt_nb(outfmt)
    organism <- .normarg_organism(organism, germline_db_name)
    germline_db_args <- .make_igblastn_germline_db_args(germline_db_name)
    germline_db_V_seqidlist <- .normarg_seqidlist(germline_db_V_seqidlist, "V")
    germline_db_D_seqidlist <- .normarg_seqidlist(germline_db_D_seqidlist, "D")
    germline_db_J_seqidlist <- .normarg_seqidlist(germline_db_J_seqidlist, "J")

    args <- c(query=query, outfmt=outfmt, organism=organism,
              germline_db_args,
              germline_db_V_seqidlist=germline_db_V_seqidlist,
              germline_db_D_seqidlist=germline_db_D_seqidlist,
              germline_db_J_seqidlist=germline_db_J_seqidlist)

    if (c_region_db_name != "") {
        c_region_db_path <- get_c_region_db_path(c_region_db_name)
        args <- c(args, c_region_db=file.path(c_region_db_path, "C"))
    }

    ## Auxiliary data.
    auxiliary_data <- file.path(igblast_root, "optional_file",
                                paste0(organism, "_gl.aux"))

    ## Turn extra args into a named character vector.
    extra_args <- list(...)
    extra_args <- setNames(as.character(extra_args), names(extra_args))

    args <- c(args, auxiliary_data=auxiliary_data, extra_args)
    cmd_args <- .make_command_line_args(args)

    if (show.command.only)
        return(.show_igblastn_command(igblast_root, cmd_args,
                                      show.in.browser=show.in.browser))

    ## Run igblastn command.
    out_file <- .run_igblastn(igblast_root, cmd_args)

    if (fmt_nb == 19) {
        ## AIRR output format is tabulated.
        AIRR_df <- read.table(out_file, header=TRUE, sep="\t")
        if (show.in.browser)
            display_data_frame_in_browser(AIRR_df)
        if (parse.out) {
            out  <- tibble(AIRR_df)
        } else {
            out <- readLines(out_file)
            class(out) <- "igblastn_raw_output"
        }
        return(out)
    }

    if (show.in.browser)
        display_local_file_in_browser(out_file)
    if (parse.out) {
        out <- .parse_igblastn_output(out_file, fmt_nb)
    } else {
        out <- readLines(out_file)
        class(out) <- "igblastn_raw_output"
    }
    out
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
    cmd_args <- if (long.help) "-help" else "-h"

    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))
    outfile <- if (show.in.browser) tempfile() else ""
    status <- system2(igblastn_exe, args=cmd_args, stdout=outfile)
    if (status != 0)
        stop(wmsg("'igblastn' returned an error"))
    if (show.in.browser)
        display_local_file_in_browser(outfile)
}


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

### Returns **absolute** path to FASTA file containing the query sequences.
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
    msg1 <- c("'outfmt' must be one of \"AIRR\", 3, 4, 7, or 19 (\"AIRR\" ",
              "is just an alias for 19), or a string describing a customized ",
              "format 7.")
    msg2 <- c("The string describing a customized format 7 must start ",
              "with a 7 followed by the desired hit table fields ",
              "(a.k.a. format specifiers) separated with spaces. ",
              "For example: \"7 std qseq sseq btop\". Note that 'std' ",
              "stands for 'qseqid sseqid pident length mismatch gapopen ",
              "gaps qstart qend sstart send evalue bitscore', which is ",
              "the default.")
    msg3 <- c("Use list_outfmt7_specifiers() to list all supported format ",
              "specifiers.")
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
    supported_specifiers <- c("std", names(list_outfmt7_specifiers()))
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
.normarg_outfmt <- function(outfmt="AIRR")
{
    if (isSingleNumber(outfmt)) {
        if (!(outfmt %in% c(3, 4, 7, 19)))
            .stop_on_invalid_outfmt()
        outfmt <- as.character(as.integer(outfmt))
    } else if (isSingleNonWhiteString(outfmt)) {
        outfmt <- trimws2(outfmt)
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

### Returns 3L, 4L, 7L, or 19L.
.extract_outfmt_nb <- function(outfmt)
{
    stopifnot(isSingleNonWhiteString(outfmt))
    outfmt_nb <- strsplit(trimws2(outfmt), " ", fixed=TRUE)[[1L]][[1L]]
    stopifnot(outfmt_nb %in% c("3", "4", "7", "19"))
    as.integer(outfmt_nb)
}

print.igblastn_raw_output <- function(x, ...) cat(x, sep="\n")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_igblastn_output()
###

### TODO: Parse output format 3 and 4.
.parse_igblastn_output <- function(out_file, outfmt_nb)
{
    stopifnot(isSingleNonWhiteString(out_file),
              isSingleInteger(outfmt_nb), outfmt_nb %in% c(3L, 4L, 7L))
    out <- readLines(out_file)
    if (outfmt_nb == 7L)
        return(parse_outfmt7(out))
    warning(wmsg("parsing of igblastn output format ",
                 outfmt_nb, " is not supported yet"))
    class(out) <- "igblastn_raw_output"
    out
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn()
###

.show_igblastn_command <- function(igblast_root, exe_args,
                                   show.in.browser=FALSE)
{
    igblastn_exe <- make_igblast_exe_path(igblast_root, "igblastn")
    cmd <- c(igblastn_exe, exe_args)
    cmd_in_1string <- paste(cmd, collapse=" ")
    outfile <- if (show.in.browser) tempfile() else ""
    cat(cmd_in_1string, "\n", file=outfile, sep="")
    if (show.in.browser)
        display_local_file_in_browser(outfile)
    cmd  # returns the command in a character vector
}

.parse_warnings_or_errors <- function(lines, pattern)
{
    stopifnot(is.character(lines), isSingleNonWhiteString(pattern))
    keep_idx <- grep(pattern, lines, ignore.case=TRUE)
    msgs <- trimws2(lines[keep_idx])
    msgs[nzchar(msgs)]
}

.parse_and_issue_warnings <- function(stderr_file)
{
    warn_msgs <- .parse_warnings_or_errors(readLines(stderr_file), "warning:")
    for (msg in warn_msgs)
        warning(wmsg(msg))
}

.stop_on_igblastn_exe_error <- function(stderr_file)
{
    err_msgs <- .parse_warnings_or_errors(readLines(stderr_file), "error:")
    if (length(err_msgs) == 0L)  # could this ever happen?
        stop(wmsg("'igblastn' returned an unknown error"))
    err_msgs <- vapply(err_msgs, wmsg, character(1), USE.NAMES=FALSE)
    stop(paste(err_msgs, collapse="\n  "))
}

### Returns the path to the output file.
.run_igblastn_exe <- function(igblast_root, exe_args)
{
    igblastn_exe <- make_igblast_exe_path(igblast_root, "igblastn")
    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))

    out_file <- tempfile()
    stderr_file <- tempfile()
    exe_args <- c(exe_args, paste("-out", out_file))
    status <- system2(igblastn_exe, args=exe_args, stderr=stderr_file)
    .parse_and_issue_warnings(stderr_file)
    if (status != 0)
        .stop_on_igblastn_exe_error(stderr_file)
    out_file
}

igblastn <- function(query, outfmt="AIRR", parse.out=TRUE, organism="auto",
                     germline_db_V="auto", germline_db_V_seqidlist=NULL,
                     germline_db_D="auto", germline_db_D_seqidlist=NULL,
                     germline_db_J="auto", germline_db_J_seqidlist=NULL,
                     c_region_db="auto", auxiliary_data=NULL, ...,
                     show.in.browser=FALSE, show.command.only=FALSE)
{
    if (!isTRUEorFALSE(parse.out))
        stop(wmsg("'parse.out' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.in.browser))
        stop(wmsg("'show.in.browser' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.command.only))
        stop(wmsg("'show.command.only' must be TRUE or FALSE"))

    igblast_root <- get_igblast_root()
    query <- .normarg_query(query)
    outfmt <- .normarg_outfmt(outfmt)
    outfmt_nb <- .extract_outfmt_nb(outfmt)
    cmd_args <- make_igblastn_command_line_args(
                     query, outfmt=outfmt, organism=organism,
                     germline_db_V=germline_db_V,
                     germline_db_V_seqidlist=germline_db_V_seqidlist,
                     germline_db_D=germline_db_D,
                     germline_db_D_seqidlist=germline_db_D_seqidlist,
                     germline_db_J=germline_db_J,
                     germline_db_J_seqidlist=germline_db_J_seqidlist,
                     c_region_db=c_region_db, auxiliary_data=auxiliary_data,
                     ...)
    exe_args <- make_exe_args(cmd_args)

    if (show.command.only) {
        out <- .show_igblastn_command(igblast_root, exe_args,
                                      show.in.browser=show.in.browser)
        return(invisible(out))
    }

    ## Run igblastn executable.
    out_file <- .run_igblastn_exe(igblast_root, exe_args)

    if (outfmt_nb == 19L) {
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
        out <- .parse_igblastn_output(out_file, outfmt_nb)
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
    exe_args <- if (long.help) "-help" else "-h"

    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))
    outfile <- if (show.in.browser) tempfile() else ""
    status <- system2(igblastn_exe, args=exe_args, stdout=outfile)
    if (status != 0)
        stop(wmsg("'igblastn' returned an error"))
    if (show.in.browser)
        display_local_file_in_browser(outfile)
}


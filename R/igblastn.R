### =========================================================================
### igblastn()
### -------------------------------------------------------------------------


### Not used at the moment.
.get_igblastr_tempdir <- function()
{
    dirpath <- file.path(tempdir(), "igblastr")
    if (!dir.exists(dirpath)) {
        if (file.exists(dirpath))
            stop(wmsg("Anomaly: '", dirpath, "' already exists ",
                      "as a file, not as a directory."))
        dir.create(dirpath)
    }
    dirpath
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_query()
###

### Returns **absolute** path to the **uncompressed** FASTA file containing
### all the query sequences. See fasta_files_as_one_uncompressed_file() in
### R/file-utils.R for the meaning of the 'safe_to_remove' attribute.
.normarg_query <- function(query)
{
    if (is.character(query))
        return(fasta_files_as_one_uncompressed_file(query, "query"))
    if (is(query, "DNAStringSet")) {
        if (is.null(names(query)))
            stop(wmsg("DNAStringSet object 'query' must have names"))
        check_seqlens(setNames(width(query), names(query)), "query")
        path <- tempfile("igblastn_query_", fileext=".fasta")
        writeXStringSet(query, path)
        attr(path, "safe_to_remove") <- TRUE
        return(path)
    }
    stop(wmsg("'query' must be a character vector that contains ",
              "the paths to the input files (FASTA), or a named ",
              "DNAStringSet object"))
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
.parse_igblastn_output <- function(out_lines, outfmt_nb)
{
    stopifnot(is.character(out_lines),
              inherits(out_lines, "igblastn_raw_output"),
              isSingleInteger(outfmt_nb), outfmt_nb %in% c(3L, 4L, 7L))
    if (outfmt_nb == 7L)
        return(parse_outfmt7(out_lines))
    warning(wmsg("parsing of igblastn output format ",
                 outfmt_nb, " is not supported yet"))
    out_lines
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn()
###

.show_igblastn_command <- function(igblast_root, exe_args,
                                   show.in.browser=FALSE)
{
    igblastn_exe <- get_igblast_exe("igblastn", igblast_root=igblast_root)
    cmd <- c(igblastn_exe, exe_args)
    cmd_in1string <- paste(cmd, collapse=" ")
    outfile <- if (show.in.browser)
               tempfile("igblastn_command_", fileext=".txt") else ""
    cat(cmd_in1string, "\n", file=outfile, sep="")
    if (show.in.browser)
        display_local_file_in_browser(outfile)
    cmd  # returns the command in a character vector
}

.parse_errors_or_warnings <- function(stderr_lines, pattern)
{
    stopifnot(is.character(stderr_lines), isSingleNonWhiteString(pattern))
    keep_idx <- grep(pattern, stderr_lines, ignore.case=TRUE)
    msgs <- trimws2(stderr_lines[keep_idx])
    msgs[nzchar(msgs)]
}

.stop_on_igblastn_exe_errors <- function(err_msgs, crash=FALSE)
{
    header <- "The 'igblastn' executable returned the following error"
    if (length(err_msgs) > 1L)
        header <- c(header, "s")
    if (crash)
        header <- c(header, " (most likely indicating a crash)")
    header <- c(header, ":")
    err_msgs <- vapply(err_msgs, wmsg, character(1), margin=4, USE.NAMES=FALSE)
    stop(wmsg(header), paste0("\n    ", err_msgs))
}

.handle_errors_or_warnings <- function(stderr_file, status, out_is_empty)
{
    stopifnot(isSingleNonWhiteString(stderr_file),
              isSingleInteger(status),
              is.logical(out_is_empty), length(out_is_empty) == 1L)
    stderr_lines <- readLines(stderr_file)
    if (length(stderr_lines) == 0L && status == 0L)
        return()
    ## When the 'igblastn' executable crashes, we typically get an empty
    ## out file. But is an empty out file always indicating a crash? We're
    ## not sure. So we add a few more conditions before concluding that we
    ## got a crash:
    ## - In our experience an empty out file and a 0 status seem to always
    ##   be associated with a crash.
    ## - An empty out file, a non-0 status, and a single stderr line
    ##   (e.g. "Segmentation fault (core dumped)") also seem to always
    ##   be associated with a crash.
    crash <- isTRUE(out_is_empty) &&
             (status == 0L || length(stderr_lines) == 1L)
    if (crash)
        .stop_on_igblastn_exe_errors(stderr_lines, crash=TRUE)
    warn_msgs <- .parse_errors_or_warnings(stderr_lines, "warning:")
    for (msg in warn_msgs)
        warning(wmsg(msg))
    if (status != 0L) {
        err_msgs <- .parse_errors_or_warnings(stderr_lines, "error:")
        if (length(err_msgs) == 0L)  # could this ever happen?
            stop(wmsg("'igblastn' returned an unknown error"))
        .stop_on_igblastn_exe_errors(err_msgs)
    }
}

### The function calls setwd() before invoking the 'igblastn' executable so
### make sure that any file path passed thru 'exe_args' (e.g. the '-out' file
### path) is an **absolute** path.
.run_igblastn_exe <- function(igblast_root, cmd_args)
{
    igblastn_exe <- get_igblast_exe("igblastn", igblast_root=igblast_root)
    exe_args <- make_exe_args(cmd_args)
    stderr_file <- tempfile("igblastn_stderr_", fileext=".txt")

    oldwd <- getwd()
    setwd(igblastr_cache(LIVE_IGDATA))
    on.exit(setwd(oldwd))
    status <- system2e(igblastn_exe, args=exe_args, stderr=stderr_file)

    out_is_empty <- NA
    if (file.exists(cmd_args$out))
        out_is_empty <- file.size(cmd_args$out) == 0
    .handle_errors_or_warnings(stderr_file, status, out_is_empty)
    unlink(stderr_file)
}

igblastn <- function(query, outfmt="AIRR",
                     germline_db_V="auto", germline_db_V_seqidlist=NULL,
                     germline_db_D="auto", germline_db_D_seqidlist=NULL,
                     germline_db_J="auto", germline_db_J_seqidlist=NULL,
                     organism="auto", c_region_db="auto",
                     custom_internal_data="auto", auxiliary_data="auto",
                     domain_system=c("imgt", "kabat"), ig_seqtype="auto",
                     ...,
                     out=NULL, parse.out=TRUE,
                     show.in.browser=FALSE, show.command.only=FALSE)
{
    domain_system <- match.arg(domain_system)
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

    ## Collect arguments that will be passed to igblastn standalone executable.
    cmd_args <- prepare_igblastn_cmdline_args(
                          query, outfmt=outfmt,
                          germline_db_V=germline_db_V,
                          germline_db_V_seqidlist=germline_db_V_seqidlist,
                          germline_db_D=germline_db_D,
                          germline_db_D_seqidlist=germline_db_D_seqidlist,
                          germline_db_J=germline_db_J,
                          germline_db_J_seqidlist=germline_db_J_seqidlist,
                          organism=organism,
                          c_region_db=c_region_db,
                          custom_internal_data=custom_internal_data,
                          auxiliary_data=auxiliary_data,
                          domain_system=domain_system,
                          ig_seqtype=ig_seqtype,
                          ...,
                          out=out)

    ## Set "safe_to_remove" files for removal on exit.
    remove_idx <- which(vapply(cmd_args,
        function(arg) isTRUE(attr(arg, "safe_to_remove")),
        logical(1)))
    if (length(remove_idx) != 0L) {
        files_to_remove_on_exit <- as.character(cmd_args[remove_idx])
        on.exit(unlink(files_to_remove_on_exit))
    }

    if (show.command.only) {
        ## Put arguments in command line format.
        exe_args <- make_exe_args(cmd_args)
        ans <- .show_igblastn_command(igblast_root, exe_args,
                                      show.in.browser=show.in.browser)
        return(invisible(ans))
    }

    ## Run the 'igblastn' standalone executable included in IgBLAST.
    .run_igblastn_exe(igblast_root, cmd_args)

    if (outfmt_nb == 19L) {
        if (show.in.browser || parse.out)
            AIRR_df <- read_igblastn_AIRR_output(cmd_args$out)
        if (show.in.browser)
            display_data_frame_in_browser(AIRR_df)
        if (parse.out) {
            ans  <- tibble(AIRR_df)
        } else {
            ans <- readLines(cmd_args$out)
            class(ans) <- "igblastn_raw_output"
        }
        return(ans)
    }

    if (show.in.browser)
        display_local_file_in_browser(cmd_args$out)
    ans <- readLines(cmd_args$out)
    class(ans) <- "igblastn_raw_output"
    if (parse.out)
        ans <- .parse_igblastn_output(ans, outfmt_nb)
    ans
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
    igblastn_exe <- get_igblast_exe("igblastn", igblast_root=igblast_root)
    exe_args <- if (long.help) "-help" else "-h"

    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))
    outfile <- file.path(tempdir(), "igblastn_help.txt")
    status <- system2e(igblastn_exe, args=exe_args, stdout=outfile)
    if (status != 0)
        stop(wmsg("'igblastn' returned an error"))
    if (show.in.browser)
        display_local_file_in_browser(outfile)
    else
        writeLines(readLines(outfile))
}


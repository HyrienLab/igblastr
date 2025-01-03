### =========================================================================
### parse_fmt7()
### -------------------------------------------------------------------------
###


### S3 generic.
query_id <- function(object) UseMethod("query_id")

.wrap_info_line <- function(info_line, width)
{
    stopifnot(isSingleNonWhiteString(info_line))
    info_prefix <- "# "
    if (has_prefix(info_line, info_prefix))
        info_line <- substr(info_line, nchar(info_prefix) + 1L,
                                       nchar(info_line))
    ## Strangely, and unintuitively, 'strwrap(x, width=n)' produces slices
    ## of max width n-1, not n, so we correct for this.
    info_lines <- strwrap(info_line, width=width+1L-nchar(info_prefix))
    paste0(info_prefix, trimws(info_lines))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_query_details()
###

.parse_query_details <- function(section_lines)
{
    stopifnot(is.character(section_lines))
    class(section_lines) <- "query_details"
    section_lines
}

query_id.query_details <- function(object)
{
    sub("^# Query: ", "", object[[2L]])
}

### Accepts 'max.nchar' argument thru the ellipsis.
summary.query_details <- function(object, ...)
{
    xargs <- list(...)
    if (length(xargs) == 0L) {
        max.nchar <- NULL
    } else {
        xargnames <- names(xargs)
        if (is.null(xargnames))
            stop(wmsg("extra arguments must be named"))
        invalid <- setdiff(xargnames, "max.nchar")
        if (length(invalid) != 0L)
            stop(wmsg("invalid argument(s): ", paste(invalid, collapse=", ")))
        max.nchar <- xargs$max.nchar
    }
    q <- query_id(object)
    if (!is.null(max.nchar) && nchar(q) > max.nchar) {
        stop <- max(max.nchar - 3L, 2L)
        q <- paste0(substr(q, 1L, stop), "...")
    }
    q
}

print.query_details <- function(x, ...) cat(x, sep="\n")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_VDJ_rearrangement_summary()
###

.parse_VDJ_rearrangement_summary <- function(section_lines)
{
    stopifnot(is.character(section_lines),
              length(section_lines) == 2L,
              identical(has_prefix(section_lines, "#"), c(TRUE, FALSE)))
    summary <- as.list(strsplit(section_lines[[2L]], "\t", fixed=TRUE)[[1L]])
    stopifnot(length(summary) == 10L)
    fields <- c("Top_V_gene_match", "Top_D_gene_match",
                "Top_J_gene_match", "Top_C_gene_match",
                "Chain_type", "stop_codon", "V_J_frame",
                "Productive", "Strand", "V_Frame_shift")
    names(summary) <- fields
    attr(summary, "info_line") <- trimws(section_lines[[1L]])
    class(summary) <- "VDJ_rearrangement_summary"
    summary
}

print.VDJ_rearrangement_summary <- function(x, ...)
{
    info_line <- attr(x, "info_line")
    lines <- .wrap_info_line(info_line, getOption("width"))
    cat(lines, sep="\n")
    cat("\n")
    attr(x, "info_line") <- NULL
    print(unclass(x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_VDJ_junction_details()
###

.parse_VDJ_junction_details <- function(section_lines)
{
    stopifnot(is.character(section_lines),
              length(section_lines) == 2L,
              identical(has_prefix(section_lines, "#"), c(TRUE, FALSE)))
    details <- as.list(strsplit(section_lines[[2L]], "\t", fixed=TRUE)[[1L]])
    stopifnot(length(details) == 5L)
    fields <- c("V_end", "V_D_junction", "D_region", "D_J_junction", "J_start")
    names(details) <- fields
    attr(details, "info_line") <- trimws(section_lines[[1L]])
    class(details) <- "VDJ_junction_details"
    details
}

print.VDJ_junction_details <- print.VDJ_rearrangement_summary


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_subregion_sequence_details()
###

.parse_subregion_sequence_details <- function(section_lines)
{
    stopifnot(is.character(section_lines),
              length(section_lines) == 2L,
              identical(has_prefix(section_lines, "#"), c(TRUE, FALSE)))
    details <- as.list(strsplit(section_lines[[2L]], "\t", fixed=TRUE)[[1L]])
    stopifnot(length(details) == 5L)
    fields <- c("sub_region", "nucleotide_sequence", "translation",
                "start", "end")
    names(details) <- fields
    attr(details, "info_line") <- trimws(section_lines[[1L]])
    class(details) <- "subregion_sequence_details"
    details
}

print.subregion_sequence_details <- print.VDJ_rearrangement_summary


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_alignment_summary()
###

.make_alignment_summary_df <- function(section_lines)
{
    fields <- c("from", "to", "length", "matches", "mismatches",
                "gaps", "percent_identity")
    table_lines <- section_lines[!has_prefix(section_lines, prefix="#")]
    file <- tempfile()
    writeLines(table_lines, file)
    df <- read.table(file, sep="\t", row.names=1L)
    colnames(df) <- fields
    df
}

.parse_alignment_summary <- function(section_lines)
{
    stopifnot(is.character(section_lines),
              which(has_prefix(section_lines, "#")) == 1L)
    summary <-.make_alignment_summary_df(section_lines)
    attr(summary, "info_line") <- trimws(section_lines[[1L]])
    class(summary) <- c("alignment_summary", class(summary))
    summary
}

print.alignment_summary <- function(x, ...)
{
    info_line <- attr(x, "info_line")
    lines <- .wrap_info_line(info_line, getOption("width"))
    cat(lines, sep="\n")
    cat("\n")
    attr(x, "info_line") <- NULL
    class(x) <- tail(class(x), n=-1L)
    print(x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_hit_table()
###

.parse_hit_table_fields <- function(section_lines)
{
    fields_line_prefix <- "# Fields: "
    fields_line <- section_lines[has_prefix(section_lines, fields_line_prefix)]
    stopifnot(length(fields_line) == 1L)
    fields_line <- substr(fields_line, nchar(fields_line_prefix) + 1L,
                                       nchar(fields_line))
    fields <- strsplit(fields_line, ", ", fixed=TRUE)[[1L]]
    ## Sanitize fields.
    fields <- chartr(" ", "_", fields)
    fields <- gsub("%", "percent", fields)
    gsub(".", "", fields, fixed=TRUE)
}

.make_hit_table_df <- function(section_lines)
{
    fields <- c("chain_type", .parse_hit_table_fields(section_lines))
    table_lines <- section_lines[!has_prefix(section_lines, prefix="#")]
    file <- tempfile()
    writeLines(table_lines, file)
    read.table(file, sep="\t", col.names=fields, check.names=FALSE)
}

.parse_hit_table <- function(section_lines)
{
    stopifnot(is.character(section_lines),
              identical(which(has_prefix(section_lines, "#")), 1:3))
    hit_table <- .make_hit_table_df(section_lines)
    attr(hit_table, "info_lines") <- trimws(section_lines[1:3])
    class(hit_table) <- c("hit_table", class(hit_table))
    hit_table
}

print.hit_table <- function(x, ...)
{
    info_lines <- attr(x, "info_lines")
    stopifnot(is.character(info_lines))
    for (i in seq_along(info_lines)) {
        lines <- .wrap_info_line(info_lines[[i]], getOption("width"))
        cat(lines, sep="\n")
    }
    cat("\n")
    attr(x, "info_lines") <- NULL
    class(x) <- tail(class(x), n=-1L)
    print(x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_fmt7record()
###

.SECTION_SEP <- ""

### Each record is expected to contain 6 sections:
###   1. query_details
###   2. VDJ_rearrangement_summary
###   3. VDJ_junction_details
###   4. subregion_sequence_details
###   5. alignment_summary
###   6. hit_table
.parse_fmt7record <- function(record_lines)
{
    record_lines <- drop_heading_and_trailing_white_lines(record_lines)
    sep_idx <- which(record_lines == .SECTION_SEP)
    section_starts <- c(1L, sep_idx + 1L)
    section_ends <- c(sep_idx - 1L, length(record_lines))
    sections <- lapply(seq_along(section_starts),
        function(i) record_lines[(section_starts[[i]]):(section_ends[[i]])])
    stopifnot(length(sections) == 6L)
    section1 <- .parse_query_details(sections[[1L]])
    section2 <- .parse_VDJ_rearrangement_summary(sections[[2L]])
    section3 <- .parse_VDJ_junction_details(sections[[3L]])
    section4 <- .parse_subregion_sequence_details(sections[[4L]])
    section5 <- .parse_alignment_summary(sections[[5L]])
    section6 <- .parse_hit_table(sections[[6L]])
    ans <- list(
        query_details=section1,
        VDJ_rearrangement_summary=section2,
        VDJ_junction_details=section3,
        subregion_sequence_details=section4,
        alignment_summary=section5,
        hit_table=section6
    )
    class(ans) <- "fmt7record"
    ans
}

query_id.fmt7record <- function(object) query_id(object$query_details)

print.fmt7record <- function(x, ...)
{
    #cat(class(x), " object\n", sep="")
    cat(class(x), " object for query id:\n", sep="")
    cat("  ", query_id(x), "\n", sep="")
    #for (section_name in names(x)) {
    #    max.nchar <- getOption("width") - nchar(section_name) - 4L
    #    cat(" $", section_name, ": ",
    #        summary(x[[section_name]], max.nchar=max.nchar), "\n", sep="")
    #}
    cat("Components:\n")
    for (section_name in names(x))
        cat("  $", section_name, "\n", sep="")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_fmt7footer()
###

.parse_fmt7footer <- function(footer_lines)
{
    stopifnot(is.character(footer_lines))
    class(footer_lines) <- "fmt7footer"
    footer_lines
}

print.fmt7footer <- function(x, ...) cat(class(x), " object\n", sep="")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### parse_fmt7()
###

.RECORD_SEP <- "# IGBLASTN"

### Returns a list with 2 components: 'records' and 'footer'.
parse_fmt7 <- function(out)
{
    footer_start <- grep("^Total queries = ", out)
    stopifnot(length(footer_start) == 1L)
    footer_lines <- out[footer_start:length(out)]
    footer <- .parse_fmt7footer(footer_lines)

    record_lines <- out[seq_len(footer_start - 1L)]
    record_starts <- which(record_lines == .RECORD_SEP)
    record_ends <- c(record_starts[-1L] - 1L, length(record_lines))
    record_ranges <- PartitioningByEnd(record_ends)
    list_of_records <- lapply(record_ranges,
        function(idx) .parse_fmt7record(record_lines[idx]))
    #class(list_of_records) <- "list_of_fmt7records"

    ans <- list(records=list_of_records, footer=footer)
    #class(ans) <- "fmt7"
    ans
}

#print.fmt7 <- function(x, ...)
#{
#    cat("igblastn ", class(x), " object with ",
#        length(x$records), " record(s) and a footer\n", sep="")
#    cat("# - access the records with <object>$records\n")
#    cat("# - access the footer with <object>$footer\n")
#    cat("# where <object> is the name of the variable containing the object\n")
#}


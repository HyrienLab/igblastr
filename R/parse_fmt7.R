### =========================================================================
### parse_fmt7()
### -------------------------------------------------------------------------
###


.wrap_comment_line <- function(comment_line, width)
{
    stopifnot(isSingleNonWhiteString(comment_line))
    comment_prefix <- "# "
    if (has_prefix(comment_line, comment_prefix))
        comment_line <- substr(comment_line, nchar(comment_prefix) + 1L,
                                             nchar(comment_line))
    ## Strangely, and unintuitively, 'strwrap(x, width=n)' produces slices
    ## of max width n-1, not n, so we correct for this.
    comment_lines <- strwrap(comment_line, width=width+1L-nchar(comment_prefix))
    paste0(comment_prefix, trimws(comment_lines))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Mapping user-specified fields (a.k.a. format specifiers) to the hit table
### fields that get effectively listed in output format 7
###

.FMT7FIELDS_IN2OUT_MAP <- c(
       qseqid="query id",
          qgi="query gi",
         qacc="query acc.",
      qaccver="query acc.ver",
         qlen="query length",
       sseqid="subject id",
    sallseqid="subject ids",
          sgi="subject gi",
       sallgi="subject gis",
         sacc="subject acc.",
      saccver="subject acc.ver",
      sallacc="subject accs",
         slen="subject length",
       qstart="q. start",
         qend="q. end",
       sstart="s. start",
         send="s. end",
         qseq="query seq",
         sseq="subject seq",
       evalue="evalue",
     bitscore="bit score",
        score="score",
       length="alignment length",
       pident="% identity",
       nident="identical",
     mismatch="mismatches",
     positive="positives",
      gapopen="gap opens",
         gaps="gaps",
         ppos="% positives",
       frames="query/sbjct frames",
       qframe="query frame",
       sframe="sbjct frame",
         btop="BTOP",
       staxid="subject tax id",
     ssciname="subject sci name",
     scomname="subject com names",
   sblastname="subject blast name",
    sskingdom="subject super kingdom",
      staxids="subject tax ids",
    sscinames="subject sci names",
    scomnames="subject com names",
  sblastnames="subject blast names",
   sskingdoms="subject super kingdoms",
       stitle="subject title",
   salltitles="subject titles",
      sstrand="subject strand",
        qcovs="% query coverage per subject",
      qcovhsp="% query coverage per hsp",
       qcovus="% query coverage per uniq subject"
)

### Exported.
list_supported_format_specifiers <- function()
{
    format_specifiers <- .FMT7FIELDS_IN2OUT_MAP
    class(format_specifiers) <- "format_specifiers"
    format_specifiers
}

print.format_specifiers <- function(x, ...)
{
    for (i in seq_along(x))
        cat(sprintf("%13s: %s", names(x)[[i]], x[[i]]), sep="\n")
}

### Translate fields from effective to user-specified.
.translate_hit_table_fields_to_user_specified <- function(fields)
{
    stopifnot(is.character(fields))
    m <- match(fields, .FMT7FIELDS_IN2OUT_MAP)
    unrecognized_idx <- which(is.na(m))
    if (length(unrecognized_idx) != 0L) {
        in1string <- paste(fields[unrecognized_idx], collapse=", ")
        stop(wmsg("unrecognized hit table field(s): ", in1string))
    }
    names(.FMT7FIELDS_IN2OUT_MAP)[m]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### qseqid() S3 generic
###

qseqid <- function(object) UseMethod("qseqid")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_query_details()
###

.parse_query_details <- function(section_lines)
{
    stopifnot(is.character(section_lines))
    class(section_lines) <- "query_details"
    section_lines
}

qseqid.query_details <- function(object)
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
    q <- qseqid(object)
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

.parse_VDJ_rearrangement_summary_fields <- function(comment_line)
{
    fields <- sub("^.* for query sequence \\(", "", comment_line)
    fields <- sub("\\).*$", "", fields)
    fields <- strsplit(fields, ", ", fixed=TRUE)[[1L]]
    ## Sanitize fields.
    chartr(" -", "__", fields)
}

.parse_VDJ_rearrangement_summary <- function(section_lines)
{
    stopifnot(is.character(section_lines),
              length(section_lines) == 2L,
              identical(has_prefix(section_lines, "#"), c(TRUE, FALSE)))
    fields <- .parse_VDJ_rearrangement_summary_fields(section_lines[[1L]])
    summary <- as.list(strsplit(section_lines[[2L]], "\t", fixed=TRUE)[[1L]])
    stopifnot(length(summary) == length(fields))
    names(summary) <- fields
    attr(summary, "comment_line") <- trimws(section_lines[[1L]])
    class(summary) <- "VDJ_rearrangement_summary"
    summary
}

print.VDJ_rearrangement_summary <- function(x, ...)
{
    comment_line <- attr(x, "comment_line")
    lines <- .wrap_comment_line(comment_line, getOption("width"))
    cat(lines, sep="\n")
    cat("\n")
    attr(x, "comment_line") <- NULL
    print(unclass(x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_VDJ_junction_details()
###

.parse_VDJ_junction_details_fields <- function(comment_line)
{
    fields <- sub("^.* based on top germline gene matches \\(", "",
                  comment_line)
    fields <- sub("\\).*$", "", fields)
    fields <- strsplit(fields, ", ", fixed=TRUE)[[1L]]
    ## Sanitize fields.
    chartr(" -", "__", fields)
}

.parse_VDJ_junction_details <- function(section_lines)
{
    stopifnot(is.character(section_lines),
              length(section_lines) == 2L,
              identical(has_prefix(section_lines, "#"), c(TRUE, FALSE)))
    fields <- .parse_VDJ_junction_details_fields(section_lines[[1L]])
    details <- as.list(strsplit(section_lines[[2L]], "\t", fixed=TRUE)[[1L]])
    stopifnot(length(details) == length(fields))
    names(details) <- fields
    attr(details, "comment_line") <- trimws(section_lines[[1L]])
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
    attr(details, "comment_line") <- trimws(section_lines[[1L]])
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
    attr(summary, "comment_line") <- trimws(section_lines[[1L]])
    class(summary) <- c("alignment_summary", class(summary))
    summary
}

print.alignment_summary <- function(x, ...)
{
    comment_line <- attr(x, "comment_line")
    lines <- .wrap_comment_line(comment_line, getOption("width"))
    cat(lines, sep="\n")
    cat("\n")
    attr(x, "comment_line") <- NULL
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
    #fields <- chartr(" ", "_", fields)
    #fields <- gsub("%", "percent", fields)
    #gsub(".", "", fields, fixed=TRUE)

    ## Translate fields from effective to user-specified.
    .translate_hit_table_fields_to_user_specified(fields)
}

.make_hit_table_df <- function(section_lines)
{
    fields <- c("chaintype", .parse_hit_table_fields(section_lines))
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
    attr(hit_table, "comment_lines") <- trimws(section_lines[1:3])
    class(hit_table) <- c("hit_table", class(hit_table))
    hit_table
}

print.hit_table <- function(x, ...)
{
    comment_lines <- attr(x, "comment_lines")
    stopifnot(is.character(comment_lines))
    for (i in seq_along(comment_lines)) {
        lines <- .wrap_comment_line(comment_lines[[i]], getOption("width"))
        cat(lines, sep="\n")
    }
    cat("\n")
    attr(x, "comment_lines") <- NULL
    class(x) <- tail(class(x), n=-1L)
    print(x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_fmt7record()
###

.SECTION_SEP <- ""

### An fmt7 record is expected to have 6 (typical) or 5 (rare) sections, or
### only 1 section if 0 hits found:
###   1. query_details
###   2. VDJ_rearrangement_summary
###   3. VDJ_junction_details
###   4. subregion_sequence_details (can be missing)
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
    stopifnot(length(sections) %in% c(1L, 5L, 6L))
    section1 <- sections[[1L]]
    if (length(sections) == 1L) {
        ## Should only happen when 0 hits found.
        last_line <- tail(section1, n=1L)
        stopifnot(last_line == "# 0 hits found")
        section1 <- head(section1, n=-1L)
        parsed_section1 <- .parse_query_details(section1)
        ans <- list(query_details=parsed_section1)
        class(ans) <- c("fmt7emptyrecord", "fmt7record")
        return(ans)
    }
    parsed_section1 <- .parse_query_details(section1)
    parsed_section2 <- .parse_VDJ_rearrangement_summary(sections[[2L]])
    parsed_section3 <- .parse_VDJ_junction_details(sections[[3L]])
    ans1 <- list(
        query_details=parsed_section1,
        VDJ_rearrangement_summary=parsed_section2,
        VDJ_junction_details=parsed_section3
    )
    if (length(sections) == 5L) {
        ## No 'subregion_sequence_details' section!
        parsed_section4 <- .parse_alignment_summary(sections[[4L]])
        parsed_section5 <- .parse_hit_table(sections[[5L]])
        ans2 <- list(alignment_summary=parsed_section4,
                     hit_table=parsed_section5)
    } else {
        parsed_section4 <- .parse_subregion_sequence_details(sections[[4L]])
        parsed_section5 <- .parse_alignment_summary(sections[[5L]])
        parsed_section6 <- .parse_hit_table(sections[[6L]])
        ans2 <- list(subregion_sequence_details=parsed_section4,
                     alignment_summary=parsed_section5,
                     hit_table=parsed_section6)
    }
    ans <- c(ans1, ans2)
    class(ans) <- "fmt7record"
    ans
}

qseqid.fmt7record <- function(object) qseqid(object$query_details)

print.fmt7record <- function(x, ...)
{
    cat(class(x)[[1L]], " object\n", sep="")
    cat("    qseqid(): ", qseqid(x), "\n", sep="")
    cat("  sections:\n")
    for (section_name in names(x))
        cat("    $", section_name, "\n", sep="")
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

print.fmt7footer <- function(x, ...) cat(class(x)[[1L]], " object\n", sep="")


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

    all_record_lines <- out[seq_len(footer_start - 1L)]
    all_record_starts <- which(all_record_lines == .RECORD_SEP)
    all_record_ends <- c(all_record_starts[-1L] - 1L, length(all_record_lines))
    all_record_ranges <- PartitioningByEnd(all_record_ends)
    list_of_records <- lapply(all_record_ranges,
        function(idx) .parse_fmt7record(all_record_lines[idx]))

    list(records=list_of_records, footer=footer)
}


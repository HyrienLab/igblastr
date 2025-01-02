### =========================================================================
### parse_igblastn_output()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normarg_outfmt()
###

normarg_outfmt <- function(outfmt="AIRR")
{
    if (isSingleString(outfmt) && toupper(outfmt) == "AIRR")
        return(19L)
    if (!(isSingleNumber(outfmt) && outfmt %in% c(3, 4, 7, 19)))
        stop(wmsg("'outfmt' must be \"AIRR\" or one of 3, 4, 7, 19 ",
                  "(19 means \"AIRR\")"))
    as.integer(outfmt)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### query() S3 generic
###

query <- function(object) UseMethod("query")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### parse_outfmt7()
###

.parse_query_details <- function(section)
{
    stopifnot(is.character(section))
    class(section) <- "query_details"
    section
}

### Exported!
query.query_details <- function(object)
{
    sub("^# Query: ", "", object[[2L]])
}

### Exported!
summary.query_details <- function(object, max.nchar=NA)
{
    q <- query(object)
    if (!is.na(max.nchar) && nchar(q) > max.nchar) {
        stop <- max(max.nchar - 3L, 2L)
        q <- paste0(substr(q, 1L, stop), "...")
    }
    q
}

.parse_VDJ_rearrangement_summary <- function(section)
{
    stopifnot(is.character(section))
    class(section) <- "VDJ_rearrangement_summary"
    section
}

### Exported!
summary.VDJ_rearrangement_summary <- function(object, ...)
{
    "bbb"
}

.parse_VDJ_junction_details <- function(section)
{
    stopifnot(is.character(section))
    class(section) <- "VDJ_junction_details"
    section
}

### Exported!
summary.VDJ_junction_details <- function(object, ...)
{
    "ccc"
}

.parse_subregion_sequence_details <- function(section)
{
    stopifnot(is.character(section))
    class(section) <- "subregion_sequence_details"
    section
}

### Exported!
summary.subregion_sequence_details <- function(object, ...)
{
    "ddd"
}

.parse_alignment_summary <- function(section)
{
    stopifnot(is.character(section))
    class(section) <- "alignment_summary"
    section
}

### Exported!
summary.alignment_summary <- function(object, ...)
{
    "eee"
}

.parse_hit_table <- function(section)
{
    stopifnot(is.character(section))
    class(section) <- "hit_table"
    section
}

### Exported!
summary.hit_table <- function(object, ...)
{
    "fff"
}

.SECTION_SEP <- ""

### Each record is expected to contain 6 sections:
###   1. query_details
###   2. VDJ_rearrangement_summary
###   3. VDJ_junction_details
###   4. subregion_sequence_details
###   5. alignment_summary
###   6. hit_table
.parse_fmt7record <- function(record)
{
    record <- drop_heading_and_trailing_white_lines(record)
    sep_idx <- which(record == .SECTION_SEP)
    section_starts <- c(1L, sep_idx + 1L)
    section_ends <- c(sep_idx - 1L, length(record))
    sections <- lapply(seq_along(section_starts),
        function(i) record[section_starts[[i]]:section_ends[[i]]])
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

### Exported!
query.fmt7record <- function(object) query(object$query_details)

### Exported!
print.fmt7record <- function(x, ...)
{
    cat(class(x), " object\n", sep="")
    for (section_name in names(x)) {
        max.nchar <- getOption("width") - nchar(section_name) - 4L
        cat("  ", section_name, ": ",
            summary(x[[section_name]], max.nchar), "\n", sep="")
    }
}

.RECORD_SEP <- "# IGBLASTN"

.parse_outfmt7 <- function(out)
{
    footer_start <- grep("^Total queries = ", out)
    stopifnot(length(footer_start) == 1L)
    footer <- out[footer_start:length(out)]

    records <- out[seq_len(footer_start - 1L)]
    record_starts <- which(records == .RECORD_SEP)
    record_ends <- c(record_starts[-1L] - 1L, length(records))
    record_ranges <- PartitioningByEnd(record_ends)
    ans <- lapply(record_ranges, function(idx) .parse_fmt7record(records[idx]))
    class(ans) <- "outfmt7"
    ans
}

### Exported!
print.outfmt7 <- function(x, ...)
{
    cat("igblastn ", class(x), " object with ",
        length(x), " record(s)\n", sep="")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### parse_igblastn_output()
###

parse_igblastn_output <- function(out, outfmt)
{
    outfmt <- normarg_outfmt(outfmt)
    if (outfmt == 7)
        return(.parse_outfmt7(out))
    #stop(wmsg("only output formats \"AIRR\" (19) or 7 ",
    #          "are supported at the moment "))
    out
}


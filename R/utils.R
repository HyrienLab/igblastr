### =========================================================================
### Various general purpose low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### TODO: wmsg() was replaced with this in S4Vectors >= 0.45.1 so drop
### wmsg2() and use wmsg() instead in BioC >= 3.23.
wmsg2 <- function(..., margin=2)
{
    width <- getOption("width") - margin
    paste0(strwrap(paste0(c(...), collapse=""), width=width),
           collapse=paste0("\n", strrep(" ", margin)))
}

### TODO: load_package_gracefully() was added to S4Vectors 0.47.6 so drop
### this and use S4Vectors:::load_package_gracefully() instead in BioC >= 3.23.
load_package_gracefully <- function(package, ...)
{
    if (!requireNamespace(package, quietly=TRUE))
        stop("Could not load package ", package, ". Is it installed?\n\n  ",
             wmsg("Note that ", ..., " requires the ", package, " package. ",
                  "Please install it with:"),
             "\n\n    BiocManager::install(\"", package, "\")")
}

### "\xc2\xa0" is some kind of weird white space that sometimes creeps
### in when scrapping dirty HTML documents found on the internet.
.WHITESPACES <- c(" ", "\t", "\r", "\n", "\xc2\xa0")

### Vectorized. Note that NAs do **not** get propagated. NA elements in 'x'
### produce FALSE elements in the output.
has_whitespace <- function(x)
{
    stopifnot(is.character(x))
    pattern <- paste0("[", paste(.WHITESPACES, collapse=""), "]")
    grepl(pattern, x, perl=TRUE)
}

### A simple wrapper to base::trimws() that starts by replacing **all**
### whitespaces in 'x' with regular spaces (" "), even non-leading
### and non-trailing whitespaces.
### Like 'base::trimws()', trimws2() is vectorized and propagates NAs.
trimws2 <- function(x)
{
    stopifnot(is.character(x))
    old <- paste(.WHITESPACES, collapse="")
    new <- strrep(" ", nchar(old))
    trimws(chartr(old, new, x), whitespace=" ")
}

### Vectorized. Note that NAs do **not** get propagated. NA elements in 'x'
### produce FALSE elements in the output.
is_white_str <- function(x) !nzchar(trimws2(x))

isSingleNonWhiteString <- function(x) isSingleString(x) && !is_white_str(x)

drop_heading_and_trailing_white_lines <- function(lines)
{
    stopifnot(is.character(lines))
    ok <- vapply(lines, is_white_str, logical(1), USE.NAMES=FALSE)
    nonwhite_idx <- which(!ok)
    if (length(nonwhite_idx) == 0L) {
        keep_idx <- integer(0)
    } else {
        keep_idx <- (nonwhite_idx[[1L]]):(nonwhite_idx[[length(nonwhite_idx)]])
    }
    lines[keep_idx]
}

### TODO: has_prefix() was added to S4Vectors 0.47.5 so drop this and
### use S4Vectors:::has_prefix() instead in BioC >= 3.23.
has_prefix <- function(x, prefix)
{
    stopifnot(is.character(x), isSingleString(prefix))
    substr(x, 1L, nchar(prefix)) == prefix
}

### TODO: has_suffix() was added to S4Vectors 0.47.5 so drop this and
### use S4Vectors:::has_suffix() instead in BioC >= 3.23.
has_suffix <- function(x, suffix)
{
    stopifnot(is.character(x), isSingleString(suffix))
    x_nc <- nchar(x)
    substr(x, x_nc - nchar(suffix) + 1L, x_nc) == suffix
}

### Not used at the moment.
.tabulate_strings_by_prefix <- function(x, prefixes)
{
    stopifnot(is.character(x), is.character(prefixes), length(prefixes) >= 1L)
    nc <- nchar(prefixes)
    stopifnot(all(nc == nc[[1L]]))
    x_prefixes <- substr(x, 1L, nc)
    m <- match(x_prefixes, prefixes)
    if (anyNA(m)) {
        in1string <- paste0(prefixes, collapse=", ")
        stop(wmsg("all strings in 'x' must start with one of ",
                  "the following prefixes: ", in1string))
    }
    setNames(tabulate(m, length(prefixes)), prefixes)
}

### Not used at the moment.
strslice <- function(x, width)
{
    stopifnot(isSingleString(x), isSingleNumber(width))
    chunks <- breakInChunks(nchar(x), chunksize=width)
    vapply(seq_along(chunks),
        function(i) substr(x, start(chunks)[i], end(chunks)[i]),
        character(1), USE.NAMES=FALSE)
}

check_seqlens <- function(seqlens, varname)
{
    stopifnot(is.integer(seqlens))
    seqids <- names(seqlens)
    stopifnot(!is.null(seqids))
    empty_idx <- which(seqlens == 0L)
    if (length(empty_idx) != 0L) {
        empty_seqids <- trimws2(seqids[empty_idx])
        if (!all(nzchar(empty_seqids)))
            stop(wmsg("some sequences in '", varname, "' are empty"))
        in1string <- paste(empty_seqids, collapse=", ")
        stop(wmsg("the following sequences in '", varname, "' ",
                  "are empty (showing seq ids): ", in1string))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### unsplit2()
###

### base::unsplit() assumes that the list elements in 'value' are vector-like
### objects that support subsetting by NAs. Unfortunately, most Vector
### derivatives in Bioconductor (e.g. IRanges, IRangesList, GRanges,
### GRangesList, DNAStringSet) do not support subsetting by NAs. So for
### example base::unsplit() doesn't work on an ordinary list of IRangesList
### objects.
### unsplit2() is a replacement for base::unsplit() that doesn't make such
### assumption. It works on any ordinary list of vector-like objects as long
### as the objects support c() and subsetting by an integer vector.
unsplit2 <- function(value, f)
{
    unlisted_value <- do.call(c, unname(value))
    f_len <- length(f)
    stopifnot(length(unlisted_value) == sum(lengths(value)),
              length(unlisted_value) == f_len)
    idx <- unlist(split(seq_len(f_len), f), use.names=FALSE)
    unlisted_value[S4Vectors:::reverseIntegerInjection(idx, f_len)]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_unique_allele_names()
###

### Generates a pool of unique suffixes made of lowercase letters.
.make_pool_of_suffixes <- function(min_pool_size)
{
    max_pool_size <- (length(letters)**8 - 1) / (length(letters) - 1) - 1
    if (min_pool_size > max_pool_size)
        stop(wmsg("too many duplicated allele names"))
    ans <- character(0)
    for (i in 1:7) {
        ans <- c(ans, mkAllStrings(letters, i))
        if (length(ans) >= min_pool_size)
            return(ans)
    }
    ## Should never happen because we checked for this condition earlier (see
    ## above).
    stop(wmsg("too many duplicated allele names"))
}

### Returns the suffixes in a named character vector parallel to 'x'.
.make_disambiguation_suffixes <- function(x)
{
    stopifnot(is.character(x))
    if (length(x) <= 1L) {
        suffixes <- character(length(x))
    } else {
        oo <- order(x)
        ordered_x <- x[oo]
        ir <- IRanges(1L, runLength(Rle(ordered_x)))
        pool_of_suffixes <- .make_pool_of_suffixes(max(width(ir)))
        suffixes <- extractList(pool_of_suffixes, ir)  # CharacterList
        suffixes[lengths(suffixes) == 1L] <- ""
        rev_oo <- S4Vectors:::reverseIntegerInjection(oo, length(oo))
        suffixes <- unlist(suffixes, use.names=FALSE)[rev_oo]
    }
    setNames(suffixes, names(x))
}

### Has its own unit tests in tests/testthat/test-utils.R!
### We make the allele names unique by adding a disambiguation suffix to
### them. Note that this is not guaranteed to make the allele names unique
### because adding the suffixes can create clashes with some other allele
### names that already looks like they have a disambiguation suffix.
### So we check for that and raise an error if it happens.
make_unique_allele_names <- function(allele_names, suffixes.only=FALSE)
{
    stopifnot(isTRUEorFALSE(suffixes.only))
    suffixes <- .make_disambiguation_suffixes(allele_names)
    ans <- add_suffix(allele_names, suffixes)
    if (anyDuplicated(ans))
        stop(wmsg("Our allele name disambiguation scheme doesn't work for ",
                  "your set of alleles, sorry. Please file an issue at ",
                  "https://github.com/HyrienLab/igblastr/issues"))
    if (suffixes.only) suffixes else ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### interweave()
###

### 'x' and 'y' must be vectors or vector-like objects of the same length.
### They must support c() and subsetting.
interweave <- function(x, y)
{
    N <- length(x)
    stopifnot(length(y) == N)
    c(x, y)[S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(N)]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### align_vectors_by_names()
###

align_vectors_by_names <- function(vectors)
{
    stopifnot(is.list(vectors), length(vectors) != 0L)
    all_names <- lapply(vectors,
        function(v) {
            nms <- names(v)
            if (is.null(nms))
                stop(wmsg("all vectors must be named"))
            if (anyDuplicated(nms))
                stop(wmsg("some vectors have duplicated names"))
            nms
        })
    unique_names <- unique(unlist(all_names))
    ans <- lapply(vectors,
        function(v) setNames(v[unique_names], unique_names))
    stopifnot(all(lengths(ans) == length(ans[[1L]])))
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rows_with_same_key_are_identical()
###

### Used in tests/testthat/test-load_auxdata.R!
rows_with_same_key_are_identical <- function(df, key)
{
    stopifnot(is.data.frame(df), isSingleNonWhiteString(key))
    rownames(df) <- NULL
    keys <- df[ , key]
    m <- match(keys, keys)
    df2 <- S4Vectors:::extract_data_frame_rows(df, m)
    identical(df, df2)
}

### The binary version of the above.
have_same_rows <- function(df1, df2, key)
{
    stopifnot(is.data.frame(df1), is.data.frame(df2),
              isSingleNonWhiteString(key))
    if (!identical(dim(df1), dim(df2)))
        return(FALSE)
    if (!identical(colnames(df1), colnames(df2)))
        return(FALSE)
    keys1 <- df1[ , key]
    keys2 <- df2[ , key]
    if (!setequal(keys1, keys2))
        return(FALSE)
    ## Some risk of false negative in the very atypical case where 'df1'
    ## and 'df2' both have rows with repeated keys but different content
    ## (i.e. rows_with_same_key_are_identical() returns FALSE on both),
    ## and the ordering of the rows below does not align them properly.
    ## For example:
    ##   df1 <- data.frame(ID=c(1,1), a=11:12)
    ##   df2 <- data.frame(ID=c(1,1), a=12:11)
    ##   have_same_rows(df1, df2, "ID")  # FALSE
    ## This will never happen if either 'df1' or 'df2' has no rows with
    ## repeated keys though.
    df1 <- S4Vectors:::extract_data_frame_rows(df1, order(keys1))
    df2 <- S4Vectors:::extract_data_frame_rows(df2, order(keys2))
    identical(df1, df2)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_discordant_rows()
### df_diff()
### complete_df_with_ref()
###

### Returns a 2-col data.frame with 1 row per discordant row pairs
### in 'df1' and 'df2'.
find_discordant_rows <- function(df1, df2, key)
{
    stopifnot(is.data.frame(df1), is.data.frame(df2),
              identical(colnames(df1), colnames(df2)),
              isSingleNonWhiteString(key))
    keys1 <- df1[ , key]
    keys2 <- df2[ , key]
    m <- match(keys1, keys2)
    mapped_idx <- which(!is.na(m))
    ok <- rep.int(TRUE, length(mapped_idx))
    for (colname in setdiff(colnames(df1), key)) {
        subcol1 <- df1[mapped_idx, colname]
        subcol2 <- df2[m[mapped_idx] , colname]
        ok <- ok & (subcol1 == subcol2 | is.na(subcol1) | is.na(subcol2))
    }
    rowidx1 <- mapped_idx[!ok]
    rowidx2 <- m[rowidx1]
    data.frame(rowidx1=rowidx1, rowidx2=rowidx2)
}

### 'rowlabels' can contain duplicates.
### Returns 'df' in a single string.
.df_as_string <- function(df, rowlabels)
{
    stopifnot(is.data.frame(df), is.character(rowlabels),
              nrow(df) == length(rowlabels))
    m <- as.matrix(format.data.frame(df, na.encode=FALSE))
    m <- cbind(names(rowlabels), rowlabels, m)
    m <- rbind(c(if (!is.null(names(rowlabels))) "", "", colnames(df)), m)
    for (j in seq_len(ncol(m)))
        m[ , j] <- format(m[ , j], justify="right")
    paste0(apply(m, 1L, paste, collapse=" "), "\n", collapse="")
}

### Returns a character vector with 1 element per row in 'disc_df1'
### (or 'disc_df2').
.rowpairs_as_strings <- function(disc_df1, disc_df2, label1, label2)
{
    stopifnot(is.data.frame(disc_df1), is.data.frame(disc_df2),
              identical(dim(disc_df1), dim(disc_df2)),
              identical(colnames(disc_df1), colnames(disc_df2)))
    vapply(seq_len(nrow(disc_df1)),
        function(i) {
            row1 <- disc_df1[i, ]
            row2 <- disc_df2[i, ]
            rowlabels <- c(rownames(row1), rownames(row2))
            names(rowlabels) <- paste0("    ", c(label1, label2), ":")
            .df_as_string(rbind(row1, row2), rowlabels)
        },
        character(1))
}

### Returns a character vector with 1 element per discordant row pair.
df_diff <- function(df1, df2, key, label1, label2)
{
    stopifnot(isSingleString(label1), isSingleString(label2))
    disc_rowpairs <- find_discordant_rows(df1, df2, key)
    disc_df1 <- df1[disc_rowpairs[[1L]], ]
    disc_df2 <- df2[disc_rowpairs[[2L]], ]
    .rowpairs_as_strings(disc_df1, disc_df2, label1, label2)
}

### Replaces missing values in 'df' with corresponding values in 'ref'.
### Raises an error that displays the differences if 'df' and 'ref' are
### discordant.
### Returns 'df' with possibly some of its missing values replaced with
### non-NA values taken from 'ref'.
complete_df_with_ref <- function(df, ref, key, what, df_label, ref_label)
{
    stopifnot(isSingleNonWhiteString(what))
    diff <- df_diff(df, ref, key, df_label, ref_label)
    if (length(diff) != 0L) {
        what <- c(df_label, " and ", ref_label, " ", what)
        msg1 <- c("  Disagreements between ", what, ":\n", diff)
        message(msg1)
        msg2 <- c(what, " contain discordant data (see ",
                  "above for the details)")
        stop(wmsg(msg2))
    }
    keys1 <- df[ , key]
    keys2 <- ref[ , key]
    m <- match(keys1, keys2)
    mapped_idx <- which(!is.na(m))
    ans <- df
    for (colname in setdiff(colnames(df), key)) {
        subcol1 <- df[mapped_idx, colname]
        subcol2 <- ref[m[mapped_idx] , colname]
        subcol1[is.na(subcol1)] <- subcol2[is.na(subcol1)]
        ans[mapped_idx, colname] <- subcol1
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### count_bin_hits()
###

### The "bins" are the adjacent integer intervals implicitly defined by
### integer vector 'bin_widths'. The values in 'bin_widths' must be
### non-negative integers that specify the "length" of each interval.
### Note that what we call the "length" (or "width") of an integer interval
### is simply the number of (integer) values in it. All the intervals are
### considered **adjacent** intervals, with the first interval starting on 1.
### So for example, setting 'bin_widths' to 'c(10L, 6L)' is defining the
### following 2 bins: the 1st bin is interval [1, 10], and the 2nd bin is
### interval [11, 16].
### Returns an integer vector parallel to 'bin_widths' where each element
### is the number of values in 'x' that fall in the corresponding bin.
### NAs in 'x' are ignored, as are values in 'x' that don't fall in any bin.
count_bin_hits <- function(x, bin_widths)
{
    stopifnot(is.integer(x), is.integer(bin_widths),
              !anyNA(bin_widths), all(bin_widths >= 0L))
    bin <- findInterval(x, c(0L, cumsum(bin_widths)), left.open=TRUE)
    setNames(tabulate(bin, nbins=length(bin_widths)), names(bin_widths))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous stuff
###

### A thin wrapper to untar() with more user-friendly error handling.
### 'exdir' should be the path to an existing directory that is
### preferrably empty.
untar2 <- function(tarfile, original_tarball_name, exdir=".")
{
    stopifnot(isSingleNonWhiteString(tarfile),
              isSingleNonWhiteString(original_tarball_name),
              isSingleNonWhiteString(exdir),
              dir.exists(exdir))
    code <- suppressWarnings(untar(tarfile, exdir=exdir))
    if (code != 0L)
        stop(wmsg("Anomaly: something went wrong during ",
                  "extraction of '", tarfile, "' (the local copy ",
                  "of '", original_tarball_name, "') to '", exdir, "'."))
}

### Returns the OS (e.g. Linux, Windows, or Darwin) and arch (e.g. x86_64
### or arm64) in a character vector of length 2, with names "OS" and "arch".
### Note that if the OS or arch cannot be obtained with Sys.info() then they
### get replaced with an NA.
get_OS_arch <- function()
{
    sys_info <- Sys.info()
    sysname <- sys_info[["sysname"]]
    if (!isSingleNonWhiteString(sysname))
        sysname <- NA_character_
    machine <- sys_info[["machine"]]
    if (!isSingleNonWhiteString(machine))
        machine <- NA_character_
    c(OS=sysname, arch=machine)
}

add_exe_suffix_on_Windows <- function(files, OS=get_OS_arch()[["OS"]])
{
    stopifnot(is.character(files), isSingleStringOrNA(OS))
    if (length(files) == 0L || is.na(OS) || !grepl("^win", tolower(OS)))
        return(files)
    paste0(files, ".exe")
}

### Returns a character vector with one string per key/val.
named_list_as_pretty_keyvals <- function(x, sep="; ")
{
    stopifnot(is.list(x))
    x_names <- names(x)
    stopifnot(!is.null(x_names))
    keys <- format(paste0(x_names, ": "))
    margin <- max(nchar(keys))
    x <- lapply(x, function(x) wmsg2(paste(x, collapse=sep), margin=margin))
    paste0(keys, as.character(x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### system2e(), system3(), and related utils
###

### Returns TRUE if environment variable BLAST_USAGE_REPORT is set to true,
### and FALSE in **any** other case (i.e. if BLAST_USAGE_REPORT is not set,
### or is set to false or gibberish).
get_igblastr_usage_report_from_BLAST_USAGE_REPORT <- function()
{
    usage_report <- Sys.getenv("BLAST_USAGE_REPORT")
    nzchar(usage_report) && isTRUE(as.logical(toupper(usage_report)))
}

system2e <- function(...)
{
    igblastr_usage_report <- getOption("igblastr_usage_report")
    old_BLAST_USAGE_REPORT <- Sys.getenv("BLAST_USAGE_REPORT")
    on.exit(Sys.setenv(BLAST_USAGE_REPORT=old_BLAST_USAGE_REPORT))
    if (isTRUE(igblastr_usage_report)) {
        Sys.setenv(BLAST_USAGE_REPORT="true")
    } else {
        Sys.setenv(BLAST_USAGE_REPORT="false")
    }
    system2(...)
}

### To use on the result of 'try(system2e(..., stdout=TRUE, stderr=TRUE))'.
system_command_worked <- function(out)
{
    if (inherits(out, "try-error"))
        return(FALSE)
    status <- attr(out, "status")
    is.null(status) || isTRUE(all.equal(status, 0L))
}

system_command_works <- function(command, args=character())
{
    out <- try(suppressWarnings(system2e(command, args=args,
                                         stdout=TRUE, stderr=TRUE)),
               silent=TRUE)
    system_command_worked(out)
}

system3 <- function(command, outfile, errfile, args=character())
{
    status <- system2e(command, args=args, stdout=outfile, stderr=errfile)
    if (file.exists(errfile)) {
        errmsg <- readLines(errfile)
        if (length(errmsg) != 0L)
            stop(paste(errmsg, collapse="\n"))
        unlink(errfile)
    }
    if (status != 0) {
        cmd_in1string <- paste(c(command, args), collapse=" ")
        stop(wmsg("command '", cmd_in1string, "' failed"))
    }
}

has_perl <- function() system_command_works("perl", args="-v")


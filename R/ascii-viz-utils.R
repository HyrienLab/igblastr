### =========================================================================
### Various low-level utilities to support ASCII visualization of alleles
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.JUSTIFY_VALUES <- c("center-left", "center-right", "left", "right")

make_ascii_double_arrow <- function(width, text="", justify=.JUSTIFY_VALUES,
                                    Larrow_char="<", Rarrow_char=">",
                                    color=NULL)
{
    stopifnot(is.integer(width), isSingleString(text),
              isSingleString(Larrow_char), nchar(Larrow_char) == 1L,
              isSingleString(Rarrow_char), nchar(Rarrow_char) == 1L)
    justify <- match.arg(justify)
    text <- rep.int(text, length(width))
    text[width < nchar(text) + 2L] <- ""
    ndashes <- width - 2L - nchar(text)
    if (justify == "center-left") {
        Ldashes <- ndashes %/% 2L
    } else if (justify == "center-right") {
        Ldashes <- (ndashes + 1L) %/% 2L
    } else if (justify == "left") {
        Ldashes <- 1L
    } else if (justify == "right") {
        Ldashes <- ndashes - 1L
    }
    Rdashes <- ndashes - Ldashes
    ans <- paste0(Larrow_char, strrep("-", pmax(Ldashes, 0L)), text,
                               strrep("-", pmax(Rdashes, 0L)), Rarrow_char)
    ans[width == 1L] <- "."
    ans[width == 0L] <- ""
    stopifnot(identical(nchar(ans), width))
    if (!is.null(color)) {
        span_tag <- sprintf("<span style=\"background: %s\">", color)
        ans <- paste0(span_tag, ans, "</span>")
    }
    ans
}

paint_XString <- function(x)
{
    stopifnot(is(x, "DNAString") || is(x, "AAString"))
    s <- as.character(x)
    class(s) <- c(seqtype(x), class(s))
    if (is(x, "DNAString"))
        return(Biostrings:::add_colors.DNA(s))
    Biostrings:::add_colors.AA(s)
}

### Returns a character matrix with fixed-width columns and same dimensions
### as the input DataFrame. The returned matrix has colnames on it and their
### widths match the widths of the columns.
### Some similarities with the makeNakedCharacterMatrixForDisplay() method
### for DataFrame objects defined in the S4Vectors.
DataFrame_as_formatted_matrix <- function(DF)
{
    stopifnot(is(DF, "DataFrame"))
    formatted_cols <- vapply(seq_along(DF),
        function(j) {
            colname <- colnames(DF)[[j]]
            col <- DF[[j]]
            justify <- if (is.character(col)) "left" else "right"
            format(c(colname, showAsCell(col)), justify=justify)
        }, character(nrow(DF) + 1L))
    colnames(formatted_cols) <- formatted_cols[1L, ]
    formatted_cols[-1L, , drop=FALSE]
}


### =========================================================================
### read_broken_table()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### matrix2df
###

matrix2df <- function(m, col2class)
{
    stopifnot(is.matrix(m), is.character(m),
              is.character(col2class), !is.null(names(col2class)),
              ncol(m) == length(col2class))
    cols <- lapply(setNames(seq_along(col2class), names(col2class)),
        function(i) {
            col <- m[ , i]  # character vector
            Class <- col2class[[i]]
            if (Class %in% c("integer", "numeric", "double"))
                col[col %in% "NA"] <- NA_character_
            as(col, Class)
        }
    )
    as.data.frame(cols, optional=TRUE, fix.empty.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_broken_table()
###

### 'x' must be a list of character vectors of variable length.
### Conceptually right-pads the list elements with empty strings ("")
### to make the list "constant-width" before unlisting it.
### Returns a character vector of length 'length(x) * width'.
right_pad_with_empty_strings_and_unlist <- function(x, width=NA)
{
    stopifnot(is.list(x), isSingleNumberOrNA(width))
    x_len <- length(x)
    if (x_len == 0L)
        return(character(0))
    x_lens <- lengths(x)
    max_x_lens <- max(x_lens)
    if (is.na(width)) {
        width <- max_x_lens
    } else {
        width <- as.integer(width)
        stopifnot(width >= max_x_lens)
    }
    y_lens <- width - x_lens
    x_seqalong <- seq_along(x)
    f <- rep.int(x_seqalong, y_lens)
    attributes(f) <- list(levels=as.character(x_seqalong), class="factor")
    y <- split(character(length(f)), f)
    collate_subscript <- rep(x_seqalong, each=2L)
    collate_subscript[2L * x_seqalong] <- x_seqalong + x_len
    unlist(c(x, y)[collate_subscript], recursive=FALSE, use.names=FALSE)
}

### Returns the table in a character matrix.
.read_jagged_table_as_character_matrix <- function(filepath)
{
    lines <- readLines(filepath)
    lines <- lines[nzchar(lines) & !has_prefix(lines, "#")]
    data <- strsplit(lines, split="[ \t]+")
    data <- right_pad_with_empty_strings_and_unlist(data)
    matrix(data, nrow=length(lines), byrow=TRUE)
}

### Read broken IgBLAST data files.
### IgBLAST data files (internal and auxiliary) are supposedly "tab-delimited"
### but they are broken in many ways:
###   - they use a mix of whitespaces for the field separators;
###   - each line contains a variable number of fields;
###   - some lines contain trailing whitespaces.
### So we cannot simply read them with read.table().
read_broken_table <- function(filepath, col2class)
{
    stopifnot(isSingleNonWhiteString(filepath),
              is.character(col2class), !is.null(names(col2class)))
    m <- .read_jagged_table_as_character_matrix(filepath)
    if (ncol(m) != length(col2class))
        stop(wmsg("error loading ", filepath, ": unexpected ",
                  "number of fields"))
    matrix2df(m, col2class)
}


### =========================================================================
### Manipulation of IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_igblast_auxiliary_data()
###

get_igblast_auxiliary_data <- function(organism)
{
    organism <- normalize_igblast_organism(organism)
    igblast_root <- get_igblast_root()
    dirpath <- file.path(igblast_root, "optional_file")
    auxiliary_data <- file.path(dirpath, paste0(organism, "_gl.aux"))
    if (!file.exists(auxiliary_data))
        stop(wmsg("no auxiliary data found in ", dirpath, " for ", organism))
    auxiliary_data
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_igblast_auxiliary_data()
###

### 'x' must be a list of character vectors of variable length.
### Conceptually right-pads the list elements with empty strings ("")
### to make the list "constant-width" before unlisting it.
### Returns a character vector of length 'length(x) * width'.
.right_pad_with_empty_strings_and_unlist <- function(x, width=NA)
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
.read_jagged_table_as_character_matrix <- function(file)
{
    lines <- readLines(file)
    lines <- lines[nzchar(lines) & !has_prefix(lines, "#")]
    data <- strsplit(lines, split="[ \t]+")
    data <- .right_pad_with_empty_strings_and_unlist(data)
    matrix(data, nrow=length(lines), byrow=TRUE)
}

.make_auxiliary_data_df_from_matrix <- function(m, filepath)
{
    if (ncol(m) != 5L)
        stop(wmsg("error loading ", filepath, ": unexpected number of fields"))
    data.frame(
        sseqid            =           m[ , 1L],
        coding_frame_start=as.integer(m[ , 2L]),
        chaintype         =           m[ , 3L],
        CDR3_stop         =as.integer(m[ , 4L]),
        extra_bps         =as.integer(m[ , 5L])
    )
}

### IgBLAST *.aux files are supposedly "tab-delimited" but they are
### broken in various ways:
###   - they use a mix of whitespaces for the field separators;
###   - each line contains a variable number of fields;
###   - some lines contain trailing whitespaces.
### So we cannot read them with read.table().
load_igblast_auxiliary_data <- function(organism)
{
    filepath <- get_igblast_auxiliary_data(organism)
    m <- .read_jagged_table_as_character_matrix(filepath)
    .make_auxiliary_data_df_from_matrix(m, filepath)
}


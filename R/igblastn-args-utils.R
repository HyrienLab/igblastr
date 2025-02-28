### =========================================================================
### Handling of igblastn() command-line arguments
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_organism()
###

### Maps 'db_name' to one of the "igblast organisms" (these are the organisms
### returned by list_igblast_organisms()).
.infer_igblast_organism_from_germline_db_name <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name))
    if (grepl("human", db_name, ignore.case=TRUE) ||
        grepl("Homo.sapiens", db_name, ignore.case=TRUE))
        return("human")
    if (grepl("mouse", db_name, ignore.case=TRUE) ||
        grepl("Mus.musculus", db_name, ignore.case=TRUE))
        return("mouse")
    if (grepl("rabbit", db_name, ignore.case=TRUE) ||
        grepl("Oryctolagus.cuniculus", db_name, ignore.case=TRUE))
        return("rabbit")
    if (grepl("rat", db_name, ignore.case=TRUE) ||
        grepl("Rattus.norvegicus", db_name, ignore.case=TRUE))
        return("rat")
    if (grepl("rhesus.monkey", db_name, ignore.case=TRUE) ||
        grepl("Macaca.mulatta", db_name, ignore.case=TRUE))
        return("rhesus_monkey")
    NA_character_
}

.normarg_organism <- function(organism="auto")
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    if (organism != "auto")
        return(normalize_igblast_organism(organism))
    germline_db_name <- use_germline_db()  # cannot be ""
    organism <- .infer_igblast_organism_from_germline_db_name(germline_db_name)
    if (is.na(organism))
        stop(wmsg("Don't know how to infer 'organism' from germline ",
                  "db name \"", germline_db_name, "\". Please set ",
                  "the 'organism' argument to the name of the IgBLAST ",
                  "internal data to use. ",
                  "Use list_igblast_organisms() to list all valid names."))
    organism
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_germline_db_X() and .normarg_c_region_db()
###

### The igblastn executable expects a "region db path" for any of
### its 'germline_db_[VDJ]' or 'c_region_db' argument. This is a path
### that does NOT point to an existing file but has a dirname() part
### that points to an existing directory.
.normalize_region_db_path <- function(region_db_path, what)
{
    stopifnot(isSingleNonWhiteString(region_db_path))
    if (file.exists(region_db_path))
        stop(wmsg(region_db_path, ": not the path to a ", what))
    dirpath <- dirname(region_db_path)
    if (!dir.exists(dirpath))
        stop(wmsg(dirpath, ": no such directory"))
    dirpath <- file_path_as_absolute(dirpath)
    file.path(dirpath, basename(region_db_path))
}

.normarg_germline_db_X <- function(germline_db_X, region_type)
{
    stopifnot(isSingleNonWhiteString(region_type))
    argname <- paste0("germline_db_", region_type)
    what <- paste0(region_type, "-region db")
    if (!isSingleNonWhiteString(germline_db_X))
        stop(wmsg("'", argname, "' must be \"auto\" or a single string ",
                  "containing the path to a ", what))
    if (germline_db_X != "auto")
        return(.normalize_region_db_path(germline_db_X, what))
    db_name <- use_germline_db()  # cannot be ""
    db_path <- get_germline_db_path(db_name)
    file.path(db_path, region_type)
}

.normarg_c_region_db <- function(c_region_db)
{
    if (is.null(c_region_db))
        return(NULL)
    what <- "C-region db"
    if (!isSingleNonWhiteString(c_region_db))
        stop(wmsg("'c_region_db' must be \"auto\", NULL, or a single string ",
                  "containing the path to a ", what))
    if (c_region_db != "auto")
        return(.normalize_region_db_path(c_region_db, what))
    db_name <- use_c_region_db()  # can be ""
    if (db_name == "")
        return(NULL)
    db_path <- get_c_region_db_path(db_name)
    file.path(db_path, "C")
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
    seqidlist <- trimws2(seqidlist)
    if (!all(nzchar(seqidlist)))
        stop(wmsg("character vector 'seqidlist' contains ",
                  "empty or white strings"))
    path <- tempfile()
    writeLines(seqidlist, path)
    path
}

.germline_db_is_internal <- function(db_path)
{
    stopifnot(isSingleNonWhiteString(db_path), dir.exists(db_path))
    db_path <- file_path_as_absolute(db_path)
    dirname(db_path) == get_germline_dbs_path()
}

.region_db_is_internal <- function(region_db_path)
{
    stopifnot(isSingleNonWhiteString(region_db_path))
    .germline_db_is_internal(dirname(region_db_path)) &&
        basename(region_db_path) %in% c("V", "D", "J")
}

### Assumes that 'region_db_path' is pointing to an **internally** managed
### region db.
.make_obtain_valid_seqids_Rcode <- function(region_db_path)
{
    stopifnot(isSingleNonWhiteString(region_db_path))
    db_path <- dirname(region_db_path)
    db_name <- basename(db_path)
    region_type <- basename(region_db_path)
    code <- sprintf("load_germline_db(\"%s\", region_types=\"%s\")",
                    db_name, region_type)
    sprintf("names(%s)", code)
}

.check_user_seqids <- function(user_seqids, region_db_path, region_type)
{
    fasta_file <- paste0(region_db_path, ".fasta")
    sequences <- readDNAStringSet(fasta_file)
    valid_seqids <- names(sequences)
    invalid_seqids <- setdiff(user_seqids, valid_seqids)
    if (length(invalid_seqids) != 0L) {
        what <- paste0(region_type, "-region db")
        in1string <- paste0(invalid_seqids, collapse=", ")
        msg1 <- c("Sequence id(s) not found in ", what, ": ", in1string)
        msg2a <- "Note that you can use:"
        code <- .make_obtain_valid_seqids_Rcode(region_db_path)
        msg2b <- "to obtain the list of valid germline sequence ids."
        stop(wmsg(msg1), "\n  ",
             wmsg(msg2a), "\n    ", code, "\n  ", wmsg(msg2b))
    }
}

.normarg_seqidlist <- function(seqidlist, region_db_path,
                               region_type=c("V", "D", "J"))
{
    region_type <- match.arg(region_type)
    if (is.null(seqidlist))
        return(NULL)
    if (inherits(seqidlist, "file")) {
        path <- .normarg_seqidlist_file(seqidlist)
    } else if (is.character(seqidlist)) {
        path <- .normarg_seqidlist_character(seqidlist)
    } else {
        stop(wmsg("'seqidlist' must be NULL, a 'file' object, ",
                  "or a character vector"))
    }
    user_seqids <- readLines(path)
    ## We only check the user-supplied sequence ids if the FASTA file
    ## associated with 'region_db_path' is guaranteed to exist and
    ## be in sync with the region db itself, which is only the case
    ## if 'region_db_path' points to an **internally** managed region db.
    if (.region_db_is_internal(region_db_path))
        .check_user_seqids(user_seqids, region_db_path, region_type)
    path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_auxiliary_data()
###

.normarg_auxiliary_data <- function(auxiliary_data)
{
    if (is.null(auxiliary_data))
        return(NULL)
    if (!isSingleNonWhiteString(auxiliary_data))
        stop(wmsg("'auxiliary_data' must be NULL or a single string ",
                  "that is the path to a file containing the coding frame ",
                  "start positions for the sequences in the J-region db"))
    file_path_as_absolute(auxiliary_data)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extra_args_as_named_character()
###

### Turns extra args into a named character vector with no NAs.
.extra_args_as_named_character <- function(...)
{
    xargs <- list(...)
    if (length(xargs) == 0L)
        return(setNames(character(0), character(0)))
    xargs_names <- names(xargs)
    if (is.null(xargs_names) || !all(nzchar(xargs_names)))
        stop(wmsg("extra arguments must be named"))
    if (any(has_whitespace(xargs_names)))
        stop(wmsg("argument names cannot contain whitespaces"))
    dupidx <- anyDuplicated(xargs_names)
    if (dupidx != 0L)
        stop(wmsg("argument '", xargs_names[[dupidx]], "' ",
                  "is defined more than once"))
    ## as.character() seems to ba able to handle any list, even nested ones,
    ## and to always return a character vector parallel to the input list.
    ## That is, it always produces a character vector with one element per
    ## top-level element in the input list. Also it doesn't propagate NAs as
    ## such but turns them into the "NA" string instead.
    setNames(as.character(xargs), xargs_names)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .check_igblastn_extra_args()
###

### 'xargs' must be a named character vector as returned by
### .extra_args_as_named_character() above.
.check_igblastn_extra_args <- function(xargs)
{
    stopifnot(is.character(xargs), !anyNA(xargs))
    xargs_names <- names(xargs)
    stopifnot(!is.null(xargs_names), !anyDuplicated(xargs_names))

    todo_url <- "https://github.com/hpages/igblastr/blob/devel/TODO"

    ## Check 'ig_seqtype' arg (see TODO file for the details).
    idx <- match("ig_seqtype", xargs_names)
    if (!is.na(idx)) {
        ig_seqtype <- xargs[[idx]]
        if (ig_seqtype == "TCR")
            stop(wmsg("'ig_seqtype=\"TCR\"' is not supported ",
                      "at the moment (see ", todo_url, ")"))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_igblastn_command_line_args()
###

### Returns the arguments in a named character vector where the names
### are valid igblastn command line argument names (e.g. "organism")
### and the values valid argument values (e.g. "rabbit").
make_igblastn_command_line_args <-
    function(query, outfmt="AIRR", organism="auto",
             germline_db_V="auto", germline_db_V_seqidlist=NULL,
             germline_db_D="auto", germline_db_D_seqidlist=NULL,
             germline_db_J="auto", germline_db_J_seqidlist=NULL,
             c_region_db="auto", auxiliary_data=NULL, ...)
{
    stopifnot(isSingleNonWhiteString(query),
              isSingleNonWhiteString(outfmt))

    organism <- .normarg_organism(organism)
    germline_db_V <- .normarg_germline_db_X(germline_db_V, "V")
    germline_db_V_seqidlist <- .normarg_seqidlist(germline_db_V_seqidlist,
                                                  germline_db_V, "V")
    germline_db_D <- .normarg_germline_db_X(germline_db_D, "D")
    germline_db_D_seqidlist <- .normarg_seqidlist(germline_db_D_seqidlist,
                                                  germline_db_D, "D")
    germline_db_J <- .normarg_germline_db_X(germline_db_J, "J")
    germline_db_J_seqidlist <- .normarg_seqidlist(germline_db_J_seqidlist,
                                                  germline_db_J, "J")
    c_region_db <- .normarg_c_region_db(c_region_db)
    auxiliary_data <- .normarg_auxiliary_data(auxiliary_data)

    cmd_args <- c(query=query, outfmt=outfmt, organism=organism,
                  germline_db_V=germline_db_V,
                  germline_db_V_seqidlist=germline_db_V_seqidlist,
                  germline_db_D=germline_db_D,
                  germline_db_D_seqidlist=germline_db_D_seqidlist,
                  germline_db_J=germline_db_J,
                  germline_db_J_seqidlist=germline_db_J_seqidlist,
                  c_region_db=c_region_db,
                  auxiliary_data=auxiliary_data)

    xargs <- .extra_args_as_named_character(...)
    .check_igblastn_extra_args(xargs)
    c(cmd_args, xargs)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_exe_args()
###

### 'cmd_args' must be a named character vector as returned by
### make_igblastn_command_line_args() above.
make_exe_args <- function(cmd_args)
{
    stopifnot(is.character(cmd_args))
    args_names <- names(cmd_args)
    stopifnot(!is.null(args_names))
    quoteme_idx <- grep(" ", cmd_args, fixed=TRUE)
    if (length(quoteme_idx) != 0L)
        cmd_args[quoteme_idx] <- paste0("'", cmd_args[quoteme_idx], "'")
    paste0("-", args_names, " ", cmd_args)
}


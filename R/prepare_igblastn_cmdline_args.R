### =========================================================================
### Handling of igblastn() command-line arguments
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### Whether 'db_path' points a cached germline db or not.
.germline_db_is_internal <- function(db_path)
{
    stopifnot(isSingleNonWhiteString(db_path))
    if (!dir.exists(db_path))
        return(FALSE)
    db_path <- file_path_as_absolute(db_path)
    dirname(db_path) == get_germline_dbs_home()  # NOT guaranteed to exist
}

### Whether 'bdb_path' **looks** like the path to an **internally** managed
### V/D/J "blast db", i.e. to a V/D/J "blast db" that belongs to a cached
### germline db.
.region_bdb_is_internal <- function(bdb_path, region_type)
{
    stopifnot(isSingleNonWhiteString(bdb_path),
              isSingleNonWhiteString(region_type))
    .germline_db_is_internal(dirname(bdb_path)) &&
        basename(bdb_path) == region_type
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_germline_db_X() and .normarg_c_region_db()
###

### The 'germline_db_[VDJ]' and 'c_region_db' arguments of the 'igblastn'
### standalone executable must be "database names" or, more precisely,
### paths to local "blast dbs".
### Note that the path to a local "blast db" must be supplied as something
### that **looks** like the path to a local file, except that, instead of
### pointing to an existing file, it points to a collection of files that
### makes up the "blast db". More precisely, the path to a local "blast db"
### looks like this:
###   1. Its dirname() part must point to the directory where the "blast db"
###      is located.
###   2. Its basename() part must be the name of an existing "blast db"
###      located in the dirname() directory. This "blast db" consists of
###      a collection of 10 binary files that were usually produced by
###      the 'makeblastdb' standalone executable. All these files are
###      expected to have a name of the form <bdb-name>.<suffix> where
###      <bdb-name> is the name of the "blast db" and <suffix> a 3 letter
###      suffix (see .BLASTDB_SUFFIXES in R/make_blastdbs.R for the list
###      of all known suffixes).
.normalize_bdb_path <- function(bdb_path, region_type)
{
    stopifnot(isSingleNonWhiteString(bdb_path))
    if (dir.exists(bdb_path)) {
        ## 'bdb_path' contains no "blast db" name so we default
        ## to 'region_type' (i.e. V, D, J, or C) for the "blast db" name.
        bdb_dirname <- bdb_path
        bdb_name <- region_type
    } else {
        bdb_dirname <- dirname(bdb_path)
        if (!dir.exists(bdb_dirname))
            stop(wmsg(bdb_dirname, ": no such directory"))
        bdb_name <- basename(bdb_path)
    }
    bdb_dirname <- file_path_as_absolute(bdb_dirname)
    file.path(bdb_dirname, bdb_name)
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
        return(.normalize_bdb_path(germline_db_X, region_type))
    db_name <- use_germline_db()  # cannot be ""
    db_path <- get_germline_db_path(db_name)
    file.path(db_path, region_type)
}

### Can return a NULL.
.normarg_c_region_db <- function(c_region_db)
{
    if (is.null(c_region_db))
        return(NULL)
    what <- "C-region db"
    if (!isSingleNonWhiteString(c_region_db))
        stop(wmsg("'c_region_db' must be \"auto\", NULL, or a single string ",
                  "containing the path to a ", what))
    if (c_region_db != "auto")
        return(.normalize_bdb_path(c_region_db, "C"))
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
    attr(path, "safe_to_remove") <- TRUE
    writeLines(seqidlist, path)
    path
}

### Assumes that 'bdb_path' is the path to an **internally** managed
### V/D/J "blast db", i.e. to a V/D/J "blast db" that belongs to a
### cached germline db.
.make_obtain_valid_seqids_Rcode <- function(bdb_path)
{
    stopifnot(isSingleNonWhiteString(bdb_path))
    db_name <- basename(dirname(bdb_path))
    region_type <- basename(bdb_path)
    code <- sprintf("load_germline_db(\"%s\", region_types=\"%s\")",
                    db_name, region_type)
    sprintf("names(%s)", code)
}

.check_user_seqids <- function(user_seqids, bdb_path, region_type)
{
    fasta_file <- paste0(bdb_path, ".fasta")
    sequences <- readDNAStringSet(fasta_file)
    valid_seqids <- names(sequences)
    invalid_seqids <- setdiff(user_seqids, valid_seqids)
    if (length(invalid_seqids) != 0L) {
        what <- paste0(region_type, "-region db")
        in1string <- paste0(invalid_seqids, collapse=", ")
        msg1 <- c("Sequence id(s) not found in ", what, ": ", in1string)
        msg2a <- "Note that you can use:"
        code <- .make_obtain_valid_seqids_Rcode(bdb_path)
        msg2b <- "to obtain the list of valid germline sequence ids."
        stop(wmsg(msg1), "\n  ",
             wmsg(msg2a), "\n    ", code, "\n  ", wmsg(msg2b))
    }
}

### Can return a NULL.
.normarg_seqidlist <- function(seqidlist,
                               bdb_path, region_type=VDJ_REGION_TYPES)
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
    ## associated with 'bdb_path' is guaranteed to exist and
    ## be in sync with the region db itself, which is only the case
    ## if 'bdb_path' points to an **internally** managed region db.
    if (.region_bdb_is_internal(bdb_path, region_type))
        .check_user_seqids(user_seqids, bdb_path, region_type)
    path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_custom_internal_data()
###

.check_user_supplied_internal_data <- function(intdata, germline_db_V, what)
{
    check_ndm_data_col2class(intdata, what=what)
    ## We will only perform the check below if the FASTA file associated
    ## with 'germline_db_V' is guaranteed to exist and be in sync with
    ## the V-region db itself, which is only the case if 'germline_db_V'
    ## points to an **internally** managed region db.
    if (!.region_bdb_is_internal(germline_db_V, "V"))
        return()
    V_fasta_file <- paste0(germline_db_V, ".fasta")
    V_names <- names(readDNAStringSet(V_fasta_file))
    num_missing <- sum(!(V_names %in% intdata[ , "allele_name"]))
    if (num_missing != 0L) {
        ## Seems slightly better than using use_germline_db().
        db_name <- basename(dirname(germline_db_V))
        warning(wmsg("Incomplete custom internal data: ", num_missing, " ",
                     "V allele(s) (out of ", length(V_names), " in ",
                     "germline db ", db_name, ") are not annotated in ", what),
                immediate.=TRUE)
    }
}

.get_auto_custom_internal_data <- function(germline_db_V, domain_system)
{
    if (!.region_bdb_is_internal(germline_db_V, "V"))
        return(NULL)
    db_path <- dirname(germline_db_V)
    intdata_path <- make_germline_db_intdata_path(db_path, FALSE, domain_system)
    if (!dir.exists(dirname(intdata_path)))  # db has no internal data dir
        return(NULL)
    if (!file.exists(intdata_path)) {
        warning(wmsg("internal data file ", basename(intdata_path), " ",
                     "not found in germline db ", basename(db_path), " ",
                     "--> using IgBLAST internal data from NCBI instead"),
                immediate.=TRUE)
        return(NULL)
    }
    intdata_path
}

.load_user_supplied_intdata <- function(custom_internal_data)
{
    if (!file.exists(custom_internal_data))
        stop(wmsg("custom FWR/CDR annotation file '",
                  custom_internal_data, "' not found"))
    if (dir.exists(custom_internal_data))
        stop(wmsg("'", custom_internal_data, "' is a directory, not a file"))
    read_ndm_data(custom_internal_data)
}

### Has its own tests in tests/testthat/test-prepare_igblastn_cmdline_args.R!
### Can return a NULL.
.normarg_custom_internal_data <- function(custom_internal_data="auto",
                                          germline_db_V, domain_system)
{
    if (is.null(custom_internal_data))
        return(NULL)
    if (is.data.frame(custom_internal_data)) {
        what <- "the supplied 'custom_internal_data' data.frame"
        .check_user_supplied_internal_data(custom_internal_data, germline_db_V,
                                           what)
        path <- tempfile()
        attr(path, "safe_to_remove") <- TRUE
        write_ndm_data(custom_internal_data, path)
        return(path)
    }
    if (!isSingleNonWhiteString(custom_internal_data))
        stop(wmsg("'custom_internal_data' must be \"auto\", or a data.frame ",
                  "containing custom FWR/CDR annotation, or the path to a ",
                  "file containing custom FWR/CDR annotation, or NULL"))
    if (custom_internal_data == "auto")
        return(.get_auto_custom_internal_data(germline_db_V, domain_system))
    intdata <- .load_user_supplied_intdata(custom_internal_data)
    what <- c("custom FWR/CDR annotation file '", custom_internal_data, "'")
    .check_user_supplied_internal_data(intdata, germline_db_V, what)
    file_path_as_absolute(custom_internal_data)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_organism()
###

### Is it really true that "any specified organism will be ignored" if
### the -custom_internal_data parameter is supplied, as they claim here?
### https://ncbi.github.io/igblast/cook/How-to-set-up.html (see Procedure
### to use custom FWR/CDR annotation). But then how is igblastn supposed to
### find the auxiliary data?
.normarg_organism <- function(organism="auto", custom_internal_data,
                              no_auto_germline_dbs)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    if (organism != "auto") {
        #if (!is.null(custom_internal_data))
        #    warning(wmsg("the specified 'organism' (", organism, ") is ",
        #                 "ignored when 'custom_internal_data' is specified"),
        #            immediate.=TRUE)
        return(normalize_igblast_organism(organism))
    }
    if (!is.null(custom_internal_data))
        return(NULL)
    if (no_auto_germline_dbs)
        stop(wmsg("'organism' must be specified when 'germline_db_V', ",
                  "'germline_db_D', and 'germline_db_J' are supplied ",
                  "(unless you also supply your own internal data thru ",
                  "the 'custom_internal_data' argument)"))
    db_name <- use_germline_db()  # cannot be ""
    organism <- infer_organism_shortname_from_db_name(db_name)
    if (is.na(organism))
        stop(wmsg("Don't know how to infer 'organism' from ",
                  "germline db name \"", db_name, "\". Please set ",
                  "the 'organism' argument to the name of the IgBLAST ",
                  "internal data to use. ",
                  "Use list_igblast_organisms() to list all valid names."))
    organism
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_auxiliary_data()
###

.check_user_supplied_auxiliary_data <- function(auxdata, germline_db_J, what)
{
    check_auxdata_col2class(auxdata, what=what)
    ## We will only perform the check below if the FASTA file associated
    ## with 'germline_db_J' is guaranteed to exist and be in sync with
    ## the J-region db itself, which is only the case if 'germline_db_J'
    ## points to an **internally** managed region db.
    if (!.region_bdb_is_internal(germline_db_J, "J"))
        return()
    J_fasta_file <- paste0(germline_db_J, ".fasta")
    J_names <- names(readDNAStringSet(J_fasta_file))
    num_missing <- sum(!(J_names %in% auxdata[ , "allele_name"]))
    if (num_missing != 0L) {
        ## Seems slightly better than using use_germline_db().
        db_name <- basename(dirname(germline_db_J))
        warning(wmsg("Incomplete auxiliary data: ", num_missing, " ",
                     "J allele(s) (out of ", length(J_names), " in ",
                     "germline db ", db_name, ") are not annotated in ", what),
                immediate.=TRUE)
    }
}

.get_auto_auxiliary_data <- function(germline_db_J)
{
    if (!.region_bdb_is_internal(germline_db_J, "J"))
        return(NULL)
    db_path <- dirname(germline_db_J)
    auxdata_path <- make_germline_db_auxdata_path(db_path)
    if (!dir.exists(dirname(auxdata_path)))  # db has no auxiliary data dir
        return(NULL)
    if (!file.exists(auxdata_path)) {
        warning(wmsg("auxiliary data file ", basename(auxdata_path), " ",
                     "not found in germline db ", basename(db_path), " ",
                     "--> using IgBLAST auxiliary data from NCBI instead"),
                immediate.=TRUE)
        return(NULL)
    }
    auxdata_path
}

.infer_organism_for_auxdata <- function(no_auto_germline_dbs)
{
    ## Note that the only way that the normalized 'organism' is NULL is if
    ## the normalized 'custom_internal_data' is **not** NULL. And the only
    ## way that this can happen is if either:
    ## - the user supplied their own internal data thru 'custom_internal_data';
    ## - or they didn't supply the 'germline_db_V' argument and the selected
    ##   germline db includes internal data.
    msg1 <- c("Failed to automatically figure out what auxiliary data ",
              "to use!")
    msg3 <- c("Please supply the path to the auxiliary data to use thru ",
              "the 'auxiliary_data' argument, or set 'auxiliary_data' ",
              "to NULL if you don't want to use any auxiliary data ",
              "(not recommended).")
    msg4 <- c("Note that you can use 'get_auxdata_path()' to obtain ",
              "the path to the auxiliary data shipped with IgBLAST ",
              "for a given organism.")
    msg5 <- c("See '?get_auxdata_path' for more information.")
    if (no_auto_germline_dbs) {
        stop(wmsg(msg1), "\n\n  ",
             wmsg(msg3), "\n  ", wmsg(msg4), "\n  ", wmsg(msg5))
    }
    db_name <- use_germline_db()  # cannot be ""
    organism <- infer_organism_shortname_from_db_name(db_name)
    if (is.na(organism)) {
        msg2 <- c("The germline db in use (", db_name, ") ",
                  "does not seem to correspond to an organism for ",
                  "which IgBLAST provides auxiliary data ",
                  "(use 'list_igblast_organisms()' to get the list of ",
                  "organisms that IgBLAST supports out of the box).")
        stop(wmsg(msg1), "\n\n  ", wmsg(msg2), "\n  ",
             wmsg(msg3), "\n  ", wmsg(msg4), "\n  ", wmsg(msg5))
    }
    organism
}

.load_user_supplied_auxdata <- function(auxiliary_data)
{
    if (!file.exists(auxiliary_data))
        stop(wmsg("auxiliary data file '", auxiliary_data, "' not found"))
    if (dir.exists(auxiliary_data))
        stop(wmsg("'", auxiliary_data, "' is a directory, not a file"))
    read_auxdata(auxiliary_data)
}

### Has its own tests in tests/testthat/test-prepare_igblastn_cmdline_args.R!
### Can return a NULL.
.normarg_auxiliary_data <- function(auxiliary_data="auto", organism,
                                    germline_db_J, no_auto_germline_dbs)
{
    if (is.null(auxiliary_data))
        return(NULL)
    if (is.data.frame(auxiliary_data)) {
        what <- "the supplied 'auxiliary_data' data.frame"
        .check_user_supplied_auxiliary_data(auxiliary_data, germline_db_J, what)
        path <- tempfile()
        attr(path, "safe_to_remove") <- TRUE
        write_auxdata(auxiliary_data, path)
        return(path)
    }
    if (!isSingleNonWhiteString(auxiliary_data))
        stop(wmsg("'auxiliary_data' must be \"auto\", or a data.frame ",
                  "containing auxiliary data, or the path to a file ",
                  "containing auxiliary data (.aux file), or NULL"))
    if (auxiliary_data == "auto") {
        auxiliary_data <- .get_auto_auxiliary_data(germline_db_J)
        if (!is.null(auxiliary_data))
            return(auxiliary_data)
        if (is.null(organism))
            organism <- .infer_organism_for_auxdata(no_auto_germline_dbs)
        return(get_auxdata_path(organism))
    }
    auxdata <- .load_user_supplied_auxdata(auxiliary_data)
    what <- c("auxiliary data file '", auxiliary_data, "'")
    .check_user_supplied_auxiliary_data(auxdata, germline_db_J, what)
    file_path_as_absolute(auxiliary_data)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_ig_seqtype()
###

.infer_ig_seqtype_from_db_name <- function(db_name)
{
    stopifnot(isSingleNonWhiteString(db_name))
    if (grepl(".IG", db_name, fixed=TRUE))
        return("Ig")
    if (grepl(".TR", db_name, fixed=TRUE))
        return("TCR")
    stop(wmsg("Don't know how to infer 'ig_seqtype' from ",
              "germline db name \"", db_name, "\". Please set ",
              "the 'ig_seqtype' argument to \"Ig\" or \"TCR\"."))
}

.normarg_ig_seqtype <- function(ig_seqtype="auto", no_auto_germline_dbs)
{
    err_msg <- "'ig_seqtype' must be \"auto\", \"Ig\", or \"TCR\""
    if (!isSingleNonWhiteString(ig_seqtype))
        stop(wmsg(err_msg))
    if (ig_seqtype != "auto") {
        t <- c(Ig="ig", Ig="bcr", TCR="tcr")
        m <- pmatch(tolower(ig_seqtype), t)
        if (is.na(m))
            stop(wmsg(err_msg))
        return(names(t)[[m]])
    }
    if (no_auto_germline_dbs)
        stop(wmsg("'ig_seqtype' must be specified when 'germline_db_V', ",
                  "'germline_db_D', and 'germline_db_J' are supplied"))
    db_name <- use_germline_db()  # cannot be ""
    .infer_ig_seqtype_from_db_name(db_name)
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

    if (FALSE) {
      todo_url <- "https://github.com/HyrienLab/igblastr/blob/devel/TODO"

      ## Check 'ig_seqtype' arg (see TODO file for the details).
      idx <- match("ig_seqtype", xargs_names)
      if (!is.na(idx)) {
        ig_seqtype <- xargs[[idx]]
        if (ig_seqtype == "TCR")
            stop(wmsg("'ig_seqtype=\"TCR\"' is not supported ",
                      "at the moment (see ", todo_url, ")"))
      }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_out()
###

### Must return an **absolute** path.
.normarg_out <- function(out)
{
    if (is.null(out)) {
        out <- tempfile("igblastn_out_", fileext=".txt")
        attr(out, "safe_to_remove") <- TRUE
        return(out)
    }
    if (!isSingleNonWhiteString(out))
        stop(wmsg("'out' must be NULL or a single (non-empty) string"))
    dirpath <- dirname(out)
    if (!dir.exists(dirpath))
        stop(wmsg(dirpath, ": no such directory"))
    dirpath <- file_path_as_absolute(dirpath)
    out <- file.path(dirpath, basename(out))
    if (file.exists(out))
        stop(wmsg(out, ": file exists"))
    out
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### prepare_igblastn_cmdline_args()
###

### Returns the arguments in a named list of single strings where the names
### are valid igblastn command line argument names (e.g. "organism") and the
### valid argument values (e.g. "rabbit").
### Note that the reason we return a list and not a character vector is
### because we sometimes need to put attributes on some of the strings
### e.g. sometimes we need to put the "safe_to_remove" attribute on
### the 'germline_db_[VDJ]_seqidlist' arguments.
prepare_igblastn_cmdline_args <-
    function(query, outfmt="AIRR",
             germline_db_V="auto", germline_db_V_seqidlist=NULL,
             germline_db_D="auto", germline_db_D_seqidlist=NULL,
             germline_db_J="auto", germline_db_J_seqidlist=NULL,
             organism="auto", c_region_db="auto",
             custom_internal_data="auto", auxiliary_data="auto",
             domain_system=c("imgt", "kabat"), ig_seqtype="auto",
             ...,
             out=NULL)
{
    stopifnot(isSingleNonWhiteString(query),
              isSingleNonWhiteString(outfmt))

    ## It will only make sense to get the name of the selected cached
    ## germline db with use_germline_db() if the user didn't supply all 3
    ## germline_db_X arguments, that is, if 'no_auto_germline_dbs' is FALSE.
    no_auto_germline_dbs <- !(identical(germline_db_V, "auto") ||
                              identical(germline_db_D, "auto") ||
                              identical(germline_db_J, "auto"))

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
    domain_system <- match.arg(domain_system)
    custom_internal_data <- .normarg_custom_internal_data(custom_internal_data,
                                                          germline_db_V,
                                                          domain_system)
    organism <- .normarg_organism(organism, custom_internal_data,
                                  no_auto_germline_dbs)
    auxiliary_data <- .normarg_auxiliary_data(auxiliary_data, organism,
                                              germline_db_J,
                                              no_auto_germline_dbs)
    ig_seqtype <- .normarg_ig_seqtype(ig_seqtype, no_auto_germline_dbs)

    cmd_args <- list(query=query, outfmt=outfmt,
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
                     ig_seqtype=ig_seqtype)

    xargs <- .extra_args_as_named_character(...)
    .check_igblastn_extra_args(xargs)

    out <- .normarg_out(out)

    S4Vectors:::delete_NULLs(c(cmd_args, as.list(xargs), list(out=out)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_exe_args()
###

### 'cmd_args' must be a named list of single strings character as
### returned by prepare_igblastn_cmdline_args() above.
### Returns an **unnamed** character vector parallel to 'cmd_args'.
make_exe_args <- function(cmd_args)
{
    stopifnot(is.list(cmd_args), all(lengths(cmd_args) == 1L))
    args_names <- names(cmd_args)
    stopifnot(!is.null(args_names))
    cmd_args <- as.character(cmd_args)
    ## For now we only put quotes around arguments that contain a space.
    ## Should we preventively quote all arguments, just in case?
    quoteme_idx <- grep(" ", cmd_args, fixed=TRUE)
    if (length(quoteme_idx) != 0L)
        cmd_args[quoteme_idx] <- shQuote(cmd_args[quoteme_idx])
    paste0("-", args_names, " ", cmd_args)
}


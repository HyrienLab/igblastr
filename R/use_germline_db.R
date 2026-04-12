### =========================================================================
### use_germline_db() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### use_germline_db()
###

.stop_on_no_selected_germline_db_yet <- function()
{
    msg <- c("You haven't selected the germline database to use ",
             "with igblastn() yet. Please select one with ",
             "use_germline_db(\"<db_name>\"). ",
             "See '?use_germline_db' for more information.")
    stop(wmsg(msg))
}

.get_selected_germline_db <- function(verbose=FALSE)
{
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    db_path <- get_db_in_use(germline_dbs_home, what="germline")
    if (db_path == "")
        .stop_on_no_selected_germline_db_yet()
    make_blastdbs(db_path, verbose=verbose)
    basename(db_path)
}

.how_to_suppress_use_germline_db_msg <- function(db_name)
{
    msg1 <- "To suppress this message, use:"
    msg2 <- c("suppressMessages(use_germline_db(\"", db_name, "\"))")
    c(wmsg(msg1), "\n    ", wmsg(msg2))
}

.notes_on_selecting_OGRDB_germline_db <- function(db_name)
{
    is_ogrdb_db <- has_prefix(db_name, "_AIRR.") ||
                   has_prefix(db_name, "_OGRDB.")
    if (!is_ogrdb_db)
        return()
    msg1 <- "If you use OGRDB in your research, please cite:"
    msg2 <- c("The current landscape of adaptive immune receptor ",
              "genomic and repertoire data: OGRDB and VDJbase")
    msg3 <- c("Lees, Peres, Klein, Amos et al., Nucleic Acids Research, ",
              "November 2025.")
    msg4 <- "https://doi.org/10.1093/nar/gkaf1094"
    message("  ", wmsg(msg1),
            "\n\n    ", wmsg(msg2, margin=4),
            "\n    ", wmsg(msg3, margin=4),
            "\n    ", wmsg(msg4, margin=4))
    is_src_db <- has_suffix(db_name, ".src")
    if (is_src_db) {
        ref_db_name <- sub("\\.src$", "", db_name)
        url <- "https://ogrdb.airr-community.org/germline_set/75"
        msg1 <- c("Use ", db_name, " only if you know what you are doing.")
        msg2 <- c("Note that the allele sequences in ", db_name, " come ",
                  "from the \"Source Set\" datasets provided by OGRDB. ",
                  "However, the AIRR-community/OGRDB maintainers ",
                  "recommend using the allele sequences from the ",
                  "\"Reference Set\" datasets for AIRR-seq analysis ",
                  "(see for example ", url, "), which are provided ",
                  "by ", ref_db_name, ".")
        message("\n  ", wmsg(msg1), "\n\n  ", wmsg(msg2))
    }
    message("\n  ", .how_to_suppress_use_germline_db_msg(db_name))
}

.note_on_selecting_IMGT_germline_db <- function(db_name)
{
    is_imgt_db <- has_prefix(db_name, "IMGT-")
    if (!is_imgt_db)
        return()
    message("  ", wmsg(IMGT_TERMS_OF_USE), "\n\n  ",
            .how_to_suppress_use_germline_db_msg(db_name))
}

.select_germline_db <- function(db_name, verbose=FALSE)
{
    check_germline_db_name(db_name)
    if (has_prefix(db_name, "_AIRR.")) {
        warn_if_old_AIRR_dbs_are_present()
    } else if (db_name %in% RENAMED_BUILTIN_GERMLINE_DBS) {
        warn_about_renamed_builtin_germline_dbs()
    } else {
        .notes_on_selecting_OGRDB_germline_db(db_name)
        .note_on_selecting_IMGT_germline_db(db_name)
    }

    db_path <- get_germline_db_path(db_name)
    make_blastdbs(db_path, verbose=verbose)

    ## Returns 'db_name' invisibly.
    set_db_in_use("germline", db_name, verbose=verbose)
}

use_germline_db <- function(db_name=NULL, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
    if (is.null(db_name))
        return(.get_selected_germline_db(verbose=verbose))
    .select_germline_db(db_name, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_germline_blastdbs()
###

### Not used at the moment and not exported!
clean_germline_blastdbs <- function()
{
    germline_dbs_home <- get_germline_dbs_home()  # NOT guaranteed to exist
    if (dir.exists(germline_dbs_home)) {
        all_db_names <- list_germline_dbs(names.only=TRUE)
        for (db_name in all_db_names) {
            db_path <- get_germline_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_germline_db()
###

.normarg_region_types <- function(region_types=NULL)
{
    if (is.null(region_types))
        return(VDJ_REGION_TYPES)
    if (!is.character(region_types) || anyNA(region_types))
        stop(wmsg("'region_types' must be NULL or ",
                  "a character vector with no NAs"))
    region_types <- toupper(region_types)
    if (length(region_types) == 1L) {
        region_types <- safeExplode(region_types)
    } else if (any(nchar(region_types) != 1L)) {
        stop(wmsg("'region_types' must have single-letter elements"))
    }
    if (!all(region_types %in% VDJ_REGION_TYPES))
        stop(wmsg("'region_types' can only contain letters V, D, or J"))
    region_types
}

### Returns the V, D, and/or J regions in a DNAStringSet object.
load_germline_db <- function(db_name, region_types=NULL)
{
    check_germline_db_name(db_name)
    db_path <- get_germline_db_path(db_name)
    region_types <- .normarg_region_types(region_types)
    db_fasta_files <- vapply(region_types,
        function(region_type) get_db_fasta_file(db_path, region_type),
        character(1), USE.NAMES=FALSE)
    readDNAStringSet(db_fasta_files)
}


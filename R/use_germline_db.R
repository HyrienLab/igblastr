### =========================================================================
### list_germline_dbs()
### -------------------------------------------------------------------------


get_germline_dbs_path <- function()
{   
    file.path(R_user_dir("igblastr", "cache"), "germline_dbs")
}

list_germline_dbs <- function()
{
    germline_dbs <- get_germline_dbs_path()
    if (!dir.exists(germline_dbs))
        return(character(0))
    all_db_names <- setdiff(list.files(germline_dbs), "USING")
    sort(all_db_names)
}

.stop_on_no_installed_germline_db_yet <- function()
{
    msg <- c("You don't have any installed germline database yet. ",
             "Use any of the install_*_germline_db() function (e.g. ",
             "install_VQUEST_germline_db()) to install at least one.")
    stop(wmsg(msg))
}

.get_germline_db <- function()
{
    all_db_names <- list_germline_dbs()
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    germline_dbs <- get_germline_dbs_path()
    using_path <- file.path(germline_dbs, "USING")
    if (!file.exists(using_path)) {
        msg <- c("You haven't selected any germline db to use yet. ",
                 "Please select one with use_germline_db(\"<db_name>\"). ",
                 "See '?use_germline_db' for more information.")
        stop(wmsg(msg))
    }
    db_name <- readLines(using_path)
    if (length(db_name) != 1L)
        stop(wmsg("Anomaly: file '", using_path, "' is corrupted."),
             "\n  ",
             wmsg("File should contain exactly one line. ",
                  "Try to repair with use_germline_db(\"<db_name>\"). ",
                  "See '?use_germline_db' for more information."))
    db_path <- file.path(germline_dbs, db_name)
    if (!dir.exists(db_path))
        stop(wmsg("Anomaly: file '", using_path, "' is invalid."),
             "\n  ",
             wmsg("File content ('", db_name, "') is not the name ",
                  "of an installed germline db. ",
                  "Try to repair with use_germline_db(\"<db_name>\"). ",
                  "See '?use_germline_db' for more information."))
    db_name
}

use_germline_db <- function(db_name=NULL)
{
    if (is.null(db_name))
        return(.get_germline_db())
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    all_db_names <- list_germline_dbs()
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    if (!(db_name %in% all_db_names)) {
        msg <- c("\"", db_name, "\" is not an installed germline db. ",
                 "Use list_germline_dbs() to list the germline ",
                 "databases that are currently installed on your ",
                 "machine (see '?list_germline_dbs'). ",
                 "Note that you can use any of the install_*_germline_db() ",
                 "function (e.g. install_VQUEST_germline_db()) to install ",
                 "additional germline databases.")
        stop(wmsg(msg))
    }
    germline_dbs <- get_germline_dbs_path()
    using_path <- file.path(germline_dbs, "USING")
    writeLines(db_name, using_path)
    db_name
}


### =========================================================================
### create_germline_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_loci_in_germline_fasta_dir()
###

.list_VDJ_fasta_files <- function(fasta_dir, loci_prefix)
{
    stopifnot(isSingleNonWhiteString(fasta_dir), dir.exists(fasta_dir),
              isSingleNonWhiteString(loci_prefix))
    pattern <- paste0("^", loci_prefix, ".[VDJ]\\.fasta$")
    fasta_files <- list_fasta_files(fasta_dir, pattern, full.names=FALSE)
    stopifnot(length(fasta_files) != 0L)
    fasta_files
}

### Returns a character vector of loci in canonical order.
.get_loci_from_germline_fasta_set <- function(fasta_files, loci_prefix)
{
    stopifnot(is.character(fasta_files),
              isSingleString(loci_prefix), loci_prefix %in% c("IG", "TR"))
    loci <- unique(sub("[VDJ]\\.fasta$", "", fasta_files))
    valid_loci <- if (loci_prefix == "IG") IG_LOCI else TR_LOCI
    stopifnot(all(loci %in% valid_loci))
    valid_loci[valid_loci %in% loci]  # return loci in canonical order
}

.check_fasta_set <- function(fasta_files, loci)
{
    stopifnot(is.character(fasta_files))
    loci2regiontypes <- map_loci_to_region_types(loci)
    for (locus in loci) {
        pattern <- paste0("^", locus)
        current_files <- grep(pattern, fasta_files, value=TRUE)
        expected_files <- paste0(locus, loci2regiontypes[[locus]], ".fasta")
        missing_files <- setdiff(expected_files, current_files)
        n <- length(missing_files)
        if (n != 0L) {
            verb <- if (n == 1L) " is" else "s are"
            in1string <- paste(missing_files, collapse=", ")
            warning(wmsg("the following file", verb, " missing ",
                         "for locus ", locus, ": ", in1string))
        }
        unexpected_files <- setdiff(current_files, expected_files)
        n <- length(unexpected_files)
        if (n != 0L) {
            verb <- if (n == 1L) " is" else "s are"
            in1string <- paste(unexpected_files, collapse=", ")
            warning(wmsg("the following file", verb, " usually not expected ",
                         "for locus ", locus, ": ", in1string))
        }
    }
}

list_loci_in_germline_fasta_dir <-
    function(fasta_dir, loci_prefix, check.fasta.set=FALSE)
{
    stopifnot(isTRUEorFALSE(check.fasta.set))
    fasta_files <- .list_VDJ_fasta_files(fasta_dir, loci_prefix)
    loci <- .get_loci_from_germline_fasta_set(fasta_files, loci_prefix)
    if (check.fasta.set)
        .check_fasta_set(fasta_files, loci)
    loci
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .collect_fasta_files()
###

.get_all_region_type_loci <- function(region_type, loci_prefix=c("IG", "TR"))
{
    stopifnot(isSingleNonWhiteString(region_type))
    region_type <- match.arg(region_type, VDJ_REGION_TYPES)
    stopifnot(isSingleNonWhiteString(loci_prefix))
    loci_prefix <- match.arg(loci_prefix)
    if (loci_prefix == "IG") {
        region_types_2_loci <- IG_REGION_TYPES_2_LOCI
    } else {
        region_types_2_loci <- TR_REGION_TYPES_2_LOCI
    }
    region_types_2_loci[[region_type]]
}

.get_loci_for_region_type <- function(region_type, selected_loci)
{
    loci_prefix <- extract_selected_loci_prefix(selected_loci)
    all_loci <- .get_all_region_type_loci(region_type, loci_prefix)
    intersect(all_loci, selected_loci)
}

.list_fasta_files_for_region_type <- function(fasta_dir,
                                              region_type=VDJ_REGION_TYPES)
{
    region_type <- match.arg(region_type)
    pattern <- paste0(region_type, "\\.fasta$")
    fasta_files <- list_fasta_files(fasta_dir, pattern, full.names=FALSE)
    if (length(fasta_files) == 0L)
        stop(wmsg("Anomaly: no ", region_type, " files found in ", fasta_dir))
    fasta_files
}

.collect_fasta_files <- function(fasta_dir, region_type, loci)
{
    wanted_loci <- .get_loci_for_region_type(region_type, loci)
    ## 'loci' should have gone thru .check_loci_for_missing_regions()
    ## so this is not supposed to happen. However it also went thru
    ## .get_effective_loci() which could have removed some loci from
    ## the original user selection. So yes, it's actually still possible
    ## that 'wanted_loci' will be empty!
    if (length(wanted_loci) == 0L)
        stop(wmsg("no fasta files found for region ", region_type, " for ",
                  "the selected loci"))
    wanted_files <- paste0(wanted_loci, region_type, ".fasta")
    found_files <- .list_fasta_files_for_region_type(fasta_dir, region_type)
    fasta_files <- intersect(wanted_files, found_files)
    if (length(fasta_files) == 0L)
        stop(wmsg("no fasta files found for region ", region_type, " for ",
                  "the selected loci"))
    missing_files <- setdiff(wanted_files, found_files)
    n <- length(missing_files)
    if (n != 0L) {
        verb <- if (n == 1L) " is" else "s are"
        in1string <- paste(missing_files, collapse=", ")
        warning(wmsg("the following file", verb, " missing ",
                     "for ", region_type, ": ", in1string), immediate.=TRUE)
    }
    file.path(fasta_dir, fasta_files)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_germline_db()
###

### Return TRUE if 'intdata' is "auto", and FALSE otherwise.
.infer_auto_intdata_from_intdata_arg <-
    function(intdata, argname="intdata", what="internal data")
{
    if (is.null(intdata) || is.data.frame(intdata))
        return(FALSE)
    if (!isSingleNonWhiteString(intdata))
        stop(wmsg("'", argname, "' must be NULL, \"auto\", or a data.frame ",
                  "containing custom ", what, ", or the path to a file ",
                  "containing custom ", what))
    intdata == "auto"
}

### Create the three "region dbs": one V-, one D-, and one J-region db.
.do_create_germline_db <- function(destdir, fasta_dir, loci,
                                   imgt.fasta.headers=FALSE,
                                   gapped=FALSE, intdata=NULL,
                                   fwrcdr_widths=IMGT_FWRCDR_WIDTHS,
                                   excluded_J_alleles=character(0),
                                   auxdata=NULL, ref_auxdata=NULL,
                                   disambiguate.allele.names=FALSE,
                                   verbose=FALSE)
{
    auto_intdata <- .infer_auto_intdata_from_intdata_arg(intdata)
    if (auto_intdata && !gapped)
        stop(wmsg("'intdata=\"auto\"' is only supported if the V allele ",
                  "sequences in the supplied FASTA files have gaps and ",
                  "'gapped' is set to TRUE"))

    auto_auxdata <- .infer_auto_intdata_from_intdata_arg(auxdata,
                                                         argname="auxdata",
                                                         what="auxiliary data")

    V_fasta_files <- .collect_fasta_files(fasta_dir, "V", loci)
    create_V_region_db(destdir, V_fasta_files,
                       gapped=gapped, auto.intdata=auto_intdata,
                       fwrcdr_widths=fwrcdr_widths,
                       disambiguate.allele.names=disambiguate.allele.names,
                       verbose=verbose)

    D_fasta_files <- .collect_fasta_files(fasta_dir, "D", loci)
    create_region_db(  destdir, D_fasta_files, region_type="D",
                       disambiguate.allele.names=disambiguate.allele.names,
                       verbose=verbose)

    J_fasta_files <- .collect_fasta_files(fasta_dir, "J", loci)
    create_J_region_db(destdir, J_fasta_files,
                       imgt.fasta.headers=imgt.fasta.headers,
                       excluded_alleles=excluded_J_alleles,
                       auto.auxdata=auto_auxdata, ref_auxdata=ref_auxdata,
                       disambiguate.allele.names=disambiguate.allele.names,
                       verbose=verbose)

    if (!(is.null(intdata) || auto_intdata)) {
        ## Add the user-supplied internal data to the db.
        if (isSingleNonWhiteString(intdata))
            intdata <- read_ndm_data(intdata)
        ## At this point 'intdata' is guaranteed to be a data.frame.
        ## TODO: Right now write_intdata_to_db() fails in a not-so-informative
        ## fashion if some V alleles in the db are not annotated in 'intdata'
        ## or if the latter has more rows than the V alleles in the db.
        ## It needs to be modified to:
        ## 1. Fail with an informative error message if some V alleles in
        ##    the db are not annotated in 'intdata'.
        ## 2. Accept an 'intdata' data.frame that has more rows than the
        ##    V alleles in the db, maybe with a warning.
        write_intdata_to_db(intdata, destdir, check.and.reorder=TRUE)
    }

    if (!(is.null(auxdata) || auto_auxdata)) {
        ## Add the user-supplied auxiliary data to the db.
        if (isSingleNonWhiteString(auxdata))
            auxdata <- read_auxdata(auxdata)
        ## At this point 'auxdata' is guaranteed to be a data.frame.
        ## TODO: write_auxdata_to_db() has the same problems as
        ## write_intdata_to_db(). See previous TODO above.
        write_auxdata_to_db(auxdata, destdir, check.and.reorder=TRUE)
    }
}

### A "germline db" is made of three "region dbs": one V-, one D-, and
### one J-region db. Uses create_region_db() and releated to create
### each "region db".
### Note that 'destdir' will typically be the path to a subdir of the
### GERMLINE_DBS cache compartment (see R/cache-utils.R for details about
### igblastr's cache organization). This subdir or any of its parent
### directories don't need to exist yet.
### See create_V_region_db() and create_J_region_db() in R/create_region_db.R
### for the roles of the 'gapped' and 'disambiguate.allele.names' arguments.
create_germline_db <- function(destdir, fasta_dir, loci,
                               imgt.fasta.headers=FALSE,
                               gapped=FALSE, intdata=NULL,
                               fwrcdr_widths=IMGT_FWRCDR_WIDTHS,
                               excluded_J_alleles=character(0),
                               auxdata=NULL, ref_auxdata=NULL,
                               disambiguate.allele.names=FALSE,
                               overwrite=FALSE, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir),
              isSingleNonWhiteString(fasta_dir), dir.exists(fasta_dir),
              isTRUEorFALSE(imgt.fasta.headers),
              isTRUEorFALSE(gapped),
              isTRUEorFALSE(disambiguate.allele.names),
              isTRUEorFALSE(overwrite),
              isTRUEorFALSE(verbose))

    if (dir.exists(destdir) && !overwrite)
        stop(wmsg(destdir, ": directory already exists. ",
                  "Use 'overwrite=TRUE' to overwrite it."))

    ## We first create the three region dbs in a temporary folder, and, only
    ## if successful, we replace 'destdir' with the temporary folder. Otherwise
    ## we destroy the temporary folder and raise an error. This achieves
    ## atomicity and avoids loosing the content of the existing 'destdir' in
    ## case something goes wrong.
    tmp_destdir <- tempfile("germline_db_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))
    .do_create_germline_db(tmp_destdir, fasta_dir, loci,
                           imgt.fasta.headers=imgt.fasta.headers,
                           gapped=gapped, intdata=intdata,
                           fwrcdr_widths=fwrcdr_widths,
                           excluded_J_alleles=excluded_J_alleles,
                           auxdata=auxdata, ref_auxdata=ref_auxdata,
                           disambiguate.allele.names=disambiguate.allele.names,
                           verbose=verbose)
    rename_file(tmp_destdir, destdir, replace=TRUE)
}


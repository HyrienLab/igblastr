### =========================================================================
### Low-level utilities to retrieve data from the AIRR-community/OGRDB
### download site
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


### See short introduction to OGRDB REST API at
### https://wordpress.vdjbase.org/index.php/ogrdb_news/downloading-germline-sets-from-the-command-line-or-api/

.OGRDB_API_URL <- "https://ogrdb.airr-community.org/api"

.encode_OGRDB_API_query <- function(query)
{
    stopifnot(is.character(query))
    if (length(query) != 0L)
        stopifnot(!is.null(names(query)))
    vapply(query, function(x) gsub("/", "%25252f", URLencode(x)), character(1))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_germline_set_from_OGRDB()
###

### Retrieves:
###     /germline/set/
###       {species}/{set_name}/{release_version}/{format}
###     or
###       {species}/{species_subgroup}/{set_name}/{release_version}/{format}
### as documented at https://ogrdb.airr-community.org/api/
###
### Typical usage:
###
###   fetch_germline_set_from_OGRDB("Human", set_name="IGH_VDJ")
###   fetch_germline_set_from_OGRDB("Mouse", species_subgroup="C57BL/6",
###                                 set_name="C57BL/6 IGH")
###
### Note that passing "Homo_sapiens" or "Mus musculus" works and produces
### the same results.
### Returns a character vector containing the nucleotide sequences in FASTA
### format, except when 'format' is set to "airr" or "airr_ex".
###
### About the "airr" and "airr_ex" formats
### --------------------------------------
###
### When 'format' is set to "airr" or "airr_ex", the function returns a
### named list containing a bunch of information about the germline set
### in adition to the nucleotide sequences.
### If 'res' is the results obtained with 'format="airr"' then all the
### sequences in the germline set are described in 'res$allele_descriptions'
### which is a data.frame with 1 row per sequence and dozens of columns.
### However, for some mysterious reason, this data.frame has more rows than
### the number of sequences returned when 'format' is not "airr" or "airr_ex".
### For example, 'res$allele_descriptions' has 245 rows for Human/IGH_VDJ
### but using 'format="ungapped"' returns only 236 sequences!
###
### Assuming 'alleles' is the data.frame obtained by removing the extra
### sequences from 'res$allele_descriptions', then the seq ids, species,
### species subgroups, loci, nucleotide sequences (ungapped and gapped),
### and region types, can be extracted with:
###
###   - seqids <- alleles$label
###
###   - species <- alleles$species$label
###
###   - species_subgroups <- alleles$species_subgroup
###
###   - loci <- alleles$locus
###
###   - region_types <- alleles$sequence_type
###
###   - ungapped_seqs <- alleles$coding_sequence
###     Note that for sequences of type V, this should match:
###       sapply(alleles$v_gene_delineations,
###         function(x) if (is.null(x)) NA_character_ else x$unaligned_sequence)
###
###   - gapped_seqs <- sapply(alleles$v_gene_delineations,
###         function(x) if (is.null(x)) NA_character_ else x$aligned_sequence)
###     Note that only sequences of type V can have gaps.
fetch_germline_set_from_OGRDB <-
    function(species, species_subgroup=NULL, set_name,
             release_version="published",
             format=c("ungapped", "gapped", "airr",
                      "ungapped_ex", "gapped_ex", "airr_ex"))
{
    stopifnot(isSingleNonWhiteString(species),
              is.null(species_subgroup) ||
                  isSingleNonWhiteString(species_subgroup),
              isSingleNonWhiteString(set_name),
              isSingleNonWhiteString(release_version))
    format <- match.arg(format)

    query <- c(species=species, species_subgroup=species_subgroup,
               set_name=set_name, release_version=release_version,
               format=format)
    query <- .encode_OGRDB_API_query(query)
    url <- paste(c(.OGRDB_API_URL, "germline/set", query), collapse="/")
    content <- getUrlContent(url, type="text", encoding="UTF-8")
    ### Is it possible that 'content' will be NA? See R/REST_API.R in
    ### the UCSC.utils package for more info.
    stopifnot(isSingleString(content))

    if (!(format %in% c("airr", "airr_ex"))) {
        fasta_lines <- strsplit(content, split="\n", fixed=TRUE)[[1L]]
        return(fasta_lines)  # nucleotide sequences in FASTA format
    }

    load_package_gracefully("jsonlite",
                            "setting 'format' to \"airr\" or \"airr_ex\"")
    parsed_json <- jsonlite::fromJSON(content, simplifyDataFrame=FALSE)
    ## Sanity checks.
    stopifnot(is.list(parsed_json),
              length(parsed_json) == 1L,
              identical(names(parsed_json), "GermlineSet"),
              is.list(parsed_json[[1L]]),
              length(parsed_json[[1L]]) == 1L)
    ans <- parsed_json[[1L]][[1L]]
    stopifnot(is.list(ans), !is.null(names(ans)),
              "allele_descriptions" %in% names(ans))
    ans$allele_descriptions <-
        jsonlite:::simplifyDataFrame(ans$allele_descriptions,
                                     flatten=FALSE, simplifyMatrix=TRUE)
    ans  # named list
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_ungapped_germline_set_from_OGRDB()
###

.infer_region_type_from_seqid <- function(seqid, locus)
{
    stopifnot(is.character(seqid), isSingleNonWhiteString(locus),
              all(has_prefix(seqid, locus)))
    pos <- nchar(locus) + 1L
    region_type <- substr(seqid, pos, pos)
    stopifnot(all(region_type %in% c("V", "D", "J")))
    region_type
}

.check_ungapped_seqs <- function(ungapped_seqs,
                                 species, species_subgroup=NULL, set_name,
                                 release_version="published",
                                 extended=FALSE)
{
    format <- if (extended) "airr_ex" else "airr"
    parsed_json <- fetch_germline_set_from_OGRDB(species,
                       species_subgroup=species_subgroup,
                       set_name=set_name,
                       release_version=release_version,
                       format=format)
    stop("paranoid.mode not ready yet")
}

### Produces between 1 to 3 FASTA files depending on 'set_name':
###   - for IGH_VDJ (Human):     IGHV.fasta, IGHD.fasta, IGHJ.fasta
###   - for IGKappa_VJ (Human):  IGKV.fasta, IGKJ.fasta
###   - for IGLambda_VJ (Human): IGLV.fasta, IGLJ.fasta
### Returns the number of files produced.
###
### To download all the germline sequences for Human:
### (see https://ogrdb.airr-community.org/germline_sets/Homo%20sapiens)
###
###   download_ungapped_germline_set_from_OGRDB("Human",
###                     set_name="IGH_VDJ", locus="IGH")
###   download_ungapped_germline_set_from_OGRDB("Human",
###                     set_name="IGKappa_VJ", locus="IGK")
###   download_ungapped_germline_set_from_OGRDB("Human",
###                     set_name="IGLambda_VJ", locus="IGL")
###
### --> produces a total of 7 FASTA files (full germline db).
###
### See download_mouse_germline_sets_from_OGRDB() below for how to
### conveniently download germline sets for Mouse.
download_ungapped_germline_set_from_OGRDB <-
    function(species, species_subgroup=NULL, set_name,
             locus=c("IGH", "IGK", "IGL"),
             release_version="published",
             extended=FALSE, destdir=".",
             paranoid.mode=FALSE)
{
    stopifnot(isTRUEorFALSE(extended),
              isSingleNonWhiteString(destdir),
              isTRUEorFALSE(paranoid.mode))
    locus <- match.arg(locus)

    format <- if (extended) "ungapped_ex" else "ungapped"
    fasta_lines <- fetch_germline_set_from_OGRDB(species,
                       species_subgroup=species_subgroup,
                       set_name=set_name,
                       release_version=release_version,
                       format=format)
    filepath <- tempfile(fileext=".fasta")
    writeLines(fasta_lines, filepath)
    ungapped_seqs <- readDNAStringSet(filepath)
    seq_region_types <- .infer_region_type_from_seqid(names(ungapped_seqs),
                                                      locus)

    if (paranoid.mode)
        .check_ungapped_seqs(ungapped_seqs,
                             species, species_subgroup=species_subgroup,
                             set_name=set_name,
                             release_version=release_version,
                             extended=extended)

    from <- paste0("germline set ", set_name, " (", species, ")")
    file_count <- 0L
    for (type in c("V", "D", "J")) {
        selected_seqs <- ungapped_seqs[seq_region_types == type]
        if (length(selected_seqs) == 0L)
            next
        filename <- paste0(locus, type, ".fasta")
        destfile <- file.path(destdir, filename)
        message("Write ", length(selected_seqs), " ", type, " regions ",
                "from ", from, " to ", filename, " ... ", appendLF=FALSE)
        if (file.exists(destfile))
            stop(wmsg(filename, " file already exists in ", destdir, "/"))
        writeXStringSet(selected_seqs, destfile)
        message("ok")
        file_count <- file_count + 1L
    }
    file_count
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_mouse_germline_sets_from_OGRDB()
###

### 'set_names' is expected to be a vector of Set Names from
### the Germline Sets table displayed at
###   https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
### except special Set Names "IGKJ (all strains)" and "IGLJ (all strains)".
### Assumes that all the Set Names in 'set_names' are of the
### form "<Species subgroup> <group>", where <group> is
### one of IGH, IGKV, an IGLV.
### Returns a 4-column matrix with colum names: set_name, species_subgroup,
### group, and locus.
.set_names_as_matrix <- function(set_names)
{
    if (!is.character(set_names) || anyNA(set_names))
        stop(wmsg("'set_names' must be a character vector with no NAs"))
    set_names <- trimws(set_names)
    if (anyDuplicated(set_names))
        stop(wmsg("'set_names' cannot contain duplicates"))
    FORBIDDEN_SET_NAMES <- c("IGKJ (all strains)", "IGLJ (all strains)")
    if (any(tolower(set_names) %in% tolower(FORBIDDEN_SET_NAMES))) {
        in1string <- paste0('"', FORBIDDEN_SET_NAMES, '"', collapse=" or ")
        stop(wmsg("'set_names' cannot contain ", in1string))
    }
    split_set_names <- strsplit(set_names, split=" ", fixed=TRUE)
    if (!all(lengths(split_set_names) == 2L))
        stop(wmsg("each Set Name in 'set_names' must contain exactly 1 space"))
    data <- unlist(split_set_names)
    if (is.null(data))
        data <- character(0)
    m <- matrix(data, ncol=2L, byrow=TRUE,
                dimnames=list(NULL, c("species_subgroup", "group")))
    loci <- substr(m[ , "group"], 1L, 3L)
    bad_idx <- which(!(loci %in% c("IGH", "IGK", "IGL")))
    if (length(bad_idx) != 0L) {
        in1string <- paste(set_names[bad_idx], collapse=", ")
        stop(wmsg("Bad Set Name(s): ", in1string, "."),
             "\n  ",
             wmsg("The first 3 letters of the 2nd part of the name ",
                  "is not a valid locus name."))
    }
    set_name <- matrix(set_names, dimnames=list(NULL, "set_name"))
    locus    <- matrix(loci, dimnames=list(NULL, "locus"))
    ans <- cbind(set_name, m, locus)
    species_subgroup <- unique(ans[ , "species_subgroup"])
    if (length(species_subgroup) > 1L)
        warning(wmsg("the supplied Set Names are from ",
                     "more than one Species subgroup"))
    ans
}

### See Germline Sets table displayed at
###   https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
### for all valid Set Names.
### Returns the number of files produced.
###
### To download all the germline sequences for Mouse strain A/J:
###
###   download_mouse_germline_sets_from_OGRDB(c("A/J IGKV", "A/J IGLV"))
###
### --> produces a total of 4 FASTA files (partial germline db).
###
### To download all the germline sequences for Mouse strain C57BL/6:
###
###   download_mouse_germline_sets_from_OGRDB("C57BL/6 IGH")
###
### --> produces a total of 5 FASTA files (partial germline db).
###
### To download all the germline sequences for Mouse strain C57BL/6J:
###
###   set_names <- c("C57BL/6J IGKV", "C57BL/6J IGLV")
###   download_mouse_germline_sets_from_OGRDB(set_names)
###
### --> produces a total of 4 FASTA files (partial germline db).
###
### To download all the germline sequences for Mouse strain CAST/EiJ:
###
###   set_names <- c("CAST/EiJ IGH", "CAST/EiJ IGKV", "CAST/EiJ IGLV")
###   download_mouse_germline_sets_from_OGRDB(set_names)
###
### --> produces a total of 7 FASTA files (full germline db).
download_mouse_germline_sets_from_OGRDB <-
    function(set_names,
             release_version="published",
             extended=FALSE, destdir=".",
             paranoid.mode=FALSE)
{
    stopifnot(isTRUEorFALSE(extended),
              isSingleNonWhiteString(destdir),
              isTRUEorFALSE(paranoid.mode))

    m <- .set_names_as_matrix(set_names)
    set_names         <- m[ , "set_name"]
    species_subgroups <- m[ , "species_subgroup"]
    loci              <- m[ , "locus"]
    file_count <- 0L
    for (i in seq_len(nrow(m))) {
        file_count <- file_count +
            download_ungapped_germline_set_from_OGRDB("Mouse",
                        species_subgroup=species_subgroups[[i]],
                        set_name=set_names[[i]], locus=loci[[i]],
                        release_version=release_version, extended=extended,
                        destdir=destdir, paranoid.mode=paranoid.mode)
    }
    file_count <- file_count +
        download_ungapped_germline_set_from_OGRDB("Mouse",
                    set_name="IGKJ (all strains)", locus="IGK",
                    release_version=release_version, extended=extended,
                    destdir=destdir, paranoid.mode=paranoid.mode)
    file_count <- file_count +
        download_ungapped_germline_set_from_OGRDB("Mouse",
                    set_name="IGLJ (all strains)", locus="IGL",
                    release_version=release_version, extended=extended,
                    destdir=destdir, paranoid.mode=paranoid.mode)
    file_count
}


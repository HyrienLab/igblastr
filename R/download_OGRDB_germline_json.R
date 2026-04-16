### =========================================================================
### download_OGRDB_germline_json() and related
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_OGRDB_germline_json()
###

.stop_if_duplicated_IGcomps <- function(IGcomps)
{
    m <- match(IGcomps, IGcomps)
    bad_idx <-which(m != seq_along(m))
    if (length(bad_idx) != 0L) {
        idx1 <- bad_idx[[1L]]
        set_name1 <- names(IGcomps)[[m[[idx1]]]]
        set_name2 <- names(IGcomps)[[idx1]]
        msg1 <- c("More than one germline set with the same \"IG* ",
                  "components\": \"", set_name1, "\" and \"", set_name2, "\" ",
                  "(both have the \"", IGcomps[[idx1]], "\" component ",
                  "in their name).")
        msg2 <- c("The specified germline sets must have unique ",
                  "\"IG* components\".")
        stop(wmsg(msg1), "\n  ", wmsg(msg2))
    }
}

download_OGRDB_germline_json <- function(organism, germline_sets,
                                         source_set=FALSE,
                                         destdir=".", overwrite=FALSE,
                                         recache=FALSE, ...)
{
    organism <- normalize_OGRDB_organism(organism)
    germline_sets <- normalize_OGRDB_germline_sets(germline_sets)
    set_names <- names(germline_sets)
    if (!isTRUEorFALSE(source_set))
        stop(wmsg("'source_set' must be TRUE or FALSE"))
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir)) {
        if (file.exists(destdir))
            stop(wmsg(destdir, ": not a directory"))
        stop(wmsg(destdir, ": no such directory"))
    }
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    ## We will download one JSON file per supplied germline set and will
    ## use the "IG* component" of each germline set name in 'set_names' to
    ## name the corresponding JSON file. Therefore, the "IG* components"
    ## of the supplied germline sets must be unique.
    IGcomps <- extract_IGcomps_from_OGRDB_set_names(organism, set_names)
    .stop_if_duplicated_IGcomps(IGcomps)

    local_files <- vapply(seq_along(germline_sets),
        function(i) {
            set_name <- set_names[[i]]
            set_version <- germline_sets[[i]]
            ## Download OGRDB germline set to local store if it's not
            ## already there.
            download_OGRDB_germline_set_to_OGRDB_store(organism,
                                    set_name, set_version,
                                    format="airr", source_set=source_set,
                                    recache=recache, ...)
        },
        character(1)
    )

    filenames <- paste0(IGcomps, ".json")  # guaranteed to be unique!
    names(local_files) <- filenames
    copy_files_to_dir(local_files, destdir, overwrite=overwrite)

    invisible(setNames(filenames, set_names))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_intdata_from_ogrdb_json()
###

### Returns an unnamed list with 1 list element per allele.
### Each list element is itself a named list (nested lists).
.extract_allele_descriptions <- function(json_path)
{
    if (!isSingleNonWhiteString(json_path))
        stop(wmsg("'json_path' must be a single (non-empty) string"))
    parsed_json <- fromJSON(json_path, simplifyDataFrame=FALSE)
    stopifnot(is.list(parsed_json), length(parsed_json) == 1L,
              identical(names(parsed_json), "GermlineSet"))
    germline_set <- parsed_json[[1L]]
    stopifnot(is.list(germline_set), length(germline_set) == 1L)
    allele_descriptions <- germline_set[[1L]]$allele_descriptions
    stopifnot(is.list(allele_descriptions),
              is.null(names(allele_descriptions)))
    allele_descriptions
}

.extract_IMGT_v_gene_delineation <- function(allele_description)
{
    stopifnot(is.list(allele_description), !is.null(names(allele_description)))
    sequence_type <- allele_description$sequence_type
    stopifnot(isSingleNonWhiteString(sequence_type),
              sequence_type %in% VDJC_REGION_TYPES)
    if (sequence_type != "V")
        return(NULL)
    allele_name <- allele_description$label
    stopifnot(isSingleNonWhiteString(allele_name))
    chain_type <- make_chain_type(sequence_type, allele_description$locus)
    v_gene_delineations <- allele_description$v_gene_delineations
    stopifnot(is.list(v_gene_delineations),
              is.null(names(v_gene_delineations)))
    for (v_gene_delineation in v_gene_delineations) {
        if (!identical(v_gene_delineation$delineation_scheme, "IMGT"))
            next
        stopifnot(all(V_GENE_DELINEATION_COLNAMES %in%
                      names(v_gene_delineation)))
        starts_ends <- v_gene_delineation[V_GENE_DELINEATION_COLNAMES]
        starts_ends <- setNames(as.character(starts_ends), names(starts_ends))
        return(c(allele_name=allele_name, starts_ends, chain_type=chain_type))
    }
    stop(wmsg("no IMGT v_gene_delineation information found ",
              "for allele ", allele_name))
}

### An R reimplementation of Python script makeogrannote.py included in
### IgBLAST. But with the important difference that we don't adjust the
### FR1/CDR1/FR2/CDR2/FR3 boundaries found in the JSON file by subtracting
### the number of gaps that precede them in the corresponding gapped V allele
### sequence, like they do in makeogrannote.py. The reason we don't do this
### is because this adjustment is meant to transform boundaries that are
### relative to the **gapped** V allele sequences into boundaries relative
### to the **ungapped** sequences, which is not necessary because the
### boundaries found in the JSON file are already relative to the **ungapped**
### V allele sequences. At least that's how they seem to be in 2026.
### Note that the makeogrannote.py script is from 2022 when OGRDB was still
### in its infancy, so maybe the early JSON files that they generated at the
### time had the FR1/CDR1/FR2/CDR2/FR3 boundaries relative to the **gapped**
### V allele sequences? Is this the reason why makeogrannote.py adjusts them?
### Returns a data.frame with 1 row per V allele in JSON file 'json_path'
### and the same columns as the data.frame returned by load_intdata() (see
### R/intdata-utils.R).
extract_intdata_from_ogrdb_json <- function(json_path)
{
    allele_descriptions <- .extract_allele_descriptions(json_path)
    data <- lapply(allele_descriptions, .extract_IMGT_v_gene_delineation)
    data <- S4Vectors:::delete_NULLs(data)
    col2class <- head(NDM_DATA_COL2CLASS, n=-1L)  # no coding_frame_start yet
    if (length(data) == 0L) {
        warning(wmsg("no V allele descriptions found ",
                     "in JSON file: ", json_path))
        data <- character(0)
    } else {
        stopifnot(identical(names(data[[1L]]), names(col2class)))
        data <- unlist(data, use.names=FALSE)
    }
    m <- matrix(data, ncol=length(col2class), byrow=TRUE)
    df <- matrix2df(m, col2class)
    cbind(df, coding_frame_start=integer(nrow(df)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_OGRDB_intdata()
###

### Returns the intdata in a data.frame.
.infer_intdata_from_OGRDB_gapped_V_sequences <-
    function(organism, germline_set, source_set=FALSE, recache=FALSE, ...)
{
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    on.exit(nuke_file(tmp_dir))
    download_OGRDB_germline_sequences(organism, germline_set,
                                      source_set=source_set,
                                      destdir=tmp_dir,
                                      recache=recache, ...)
    V_fasta_file <- list_fasta_files(tmp_dir, pattern="V\\.fasta$")
    if (length(V_fasta_file) == 0L) {
        msg <- c(organism, " germline set \"", names(germline_set), "\" ",
                 "(version ", germline_set, ") does not seem to contain ",
                 "any V allele")
        stop(wmsg(msg))
    }
    stopifnot(length(V_fasta_file) == 1L)

    intdata <- compute_V_gene_delineations(V_fasta_file)
    locus <- extract_loci_from_OGRDB_set_names(organism, names(germline_set))
    intdata$chain_type <- make_chain_type("V", locus)
    intdata[ , names(NDM_DATA_COL2CLASS)]
}

### Returns the intdata in a data.frame.
.fetch_intdata_from_OGRDB <-
    function(organism, germline_set, source_set=FALSE, recache=FALSE, ...)
{
    ## Download OGRDB germline set to local store if it's not already there.
    local_file <- download_OGRDB_germline_set_to_OGRDB_store(organism,
                                 names(germline_set), germline_set,
                                 format="airr", source_set=source_set,
                                 recache=recache, ...)
    extract_intdata_from_ogrdb_json(local_file)
}

### Not exported!
### Validates the intdata associated with a given OGRDB germline set by
### comparing the two methods of acquisition:
###   1. Infer intdata from the gaps in the V allele sequences.
###   2. Extract intdata from OGRDB json file.
### Returns TRUE if the two methods produce exactly the same intdata, or
### FALSE if they don't.
validate_OGRDB_intdata <- function(organism, germline_set, source_set=FALSE,
                                   recache=FALSE, ...)
{
    organism <- normalize_OGRDB_organism(organism)
    germline_set <- normalize_OGRDB_germline_sets(germline_set)
    stopifnot(length(germline_set) == 1L)
    if (!isTRUEorFALSE(source_set))
        stop(wmsg("'source_set' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    ## Method 1: Infer intdata from the gaps in the V allele sequences.
    intdata1 <- .infer_intdata_from_OGRDB_gapped_V_sequences(
                                    organism, germline_set,
                                    source_set=source_set,
                                    recache=recache, ...)


    ## Method 2: Extract intdata from OGRDB json file.
    intdata2 <- .fetch_intdata_from_OGRDB(
                                    organism, germline_set,
                                    source_set=source_set,
                                    recache=recache, ...)

    ## Compare.
    same_ndm_data(intdata1, intdata2)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_auxdata_from_ogrdb_json()
###

.extract_J_allele_length <- function(allele_description)
{
    stopifnot(is.list(allele_description), !is.null(names(allele_description)))
    coding_sequence <- allele_description$coding_sequence
    stopifnot(isSingleNonWhiteString(coding_sequence))
    seq_len <- nchar(coding_sequence)
    ## Some sanity checks.
    gene_start <- allele_description$gene_start
    gene_end <- allele_description$gene_end
    stopifnot(isSingleInteger(gene_start), isSingleInteger(gene_end),
              gene_end - gene_start + 1L == seq_len)
    seq_len
}

### IMPORTANT NOTE about the j_cdr3_end reported in OGRDB json files:
### The j_cdr3_end reported in OGRDB json files seems to be 2-based i.e.
### it's set to the 1-based value + 1 or to the 0-based value + 2.
### For example, for human allele IGHJ1*01, the json file reports:
###
###   j_codon_frame = 1 and j_cdr3_end = 19
###
### IGHJ1*01 is the first J allele annotated in human_gl.aux in IgBLAST and
### the file entry for this allele is:
###
###   allele name | first coding | chain type | CDR3 stop |  extra bps
###               |  frame start |            |           |   beyond J
###               |     position |            |           | coding end
###   ------------|--------------|------------|-----------|-----------
###   IGHJ1*01    |            0 | JH         |        17 |          1
###
### All positions in human_gl.aux are 0-based so, according to the above
### entry, the 1-based position of the CDR3 end would be 18.
### Given that j_codon_frame is 1, the 1-based CDR3 end is a multiple of 3,
### as expected.
### Furthermore, the coding sequence for IGHJ1*01 is:
###
###   GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG
###
### which translates to AEYFQHWGQGTLVTVSS if we consider that the first 3
### nucleotides of the sequence form the first codon (as suggested by the
### value of j_codon_frame).
### As expected, amino acids WGQG at position 7-10 match the WGXG motif that
### marks the start of the FWR4 on the heavy chain. This would mean that the
### first 6 amino acids (AEYFQH) are the last 6 amino acids of the CDR3,
### confirming that the 18th nucleotide in the coding sequence is where the
### CDR3 ends.
.extract_j_annotation <- function(allele_description)
{
    stopifnot(is.list(allele_description), !is.null(names(allele_description)))
    sequence_type <- allele_description$sequence_type
    stopifnot(isSingleNonWhiteString(sequence_type),
              sequence_type %in% VDJC_REGION_TYPES)
    if (sequence_type != "J")
        return(NULL)
    allele_name <- allele_description$label
    stopifnot(isSingleNonWhiteString(allele_name))
    chain_type <- make_chain_type(sequence_type, allele_description$locus)
    j_codon_frame <- allele_description$j_codon_frame
    j_cdr3_end <- allele_description$j_cdr3_end
    stopifnot(isSingleInteger(j_codon_frame), isSingleInteger(j_cdr3_end))
    coding_frame_start <- j_codon_frame - 1L
    seq_len <- .extract_J_allele_length(allele_description)
    extra_bps <- (seq_len - coding_frame_start) %% 3L
    c(allele_name=allele_name,
      coding_frame_start=coding_frame_start,
      chain_type=chain_type,
      cdr3_end=j_cdr3_end - 2L,  # see IMPORTANT NOTE above
      extra_bps=extra_bps)
}

extract_auxdata_from_ogrdb_json <- function(json_path)
{
    allele_descriptions <- .extract_allele_descriptions(json_path)
    data <- lapply(allele_descriptions, .extract_j_annotation)
    data <- S4Vectors:::delete_NULLs(data)
    if (length(data) == 0L) {
        warning(wmsg("no J allele descriptions found ",
                     "in JSON file: ", json_path))
        data <- character(0)
    } else {
        stopifnot(identical(names(data[[1L]]), names(AUXDATA_COL2CLASS)))
        data <- unlist(data, use.names=FALSE)
    }
    m <- matrix(data, ncol=length(AUXDATA_COL2CLASS), byrow=TRUE)
    matrix2df(m, AUXDATA_COL2CLASS)
}


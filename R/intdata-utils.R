### =========================================================================
### Access and manipulate IgBLAST internal data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_germline_db_intdata_path()
### get_intdata_path()
###

.make_intdata_file_suffix <- function(for.aa, domain_system)
{
    stopifnot(isTRUEorFALSE(for.aa), isSingleNonWhiteString(domain_system))
    paste0(".", if (for.aa) "pdm" else "ndm", ".", domain_system)
}

.get_igblast_intdata_path <- function(which, organism, for.aa, domain_system)
{
    organism <- normalize_igblast_organism(organism)
    intdata_dir <- file.path(path_to_igdata(which), "internal_data", organism)
    if (!dir.exists(intdata_dir))
        stop(wmsg("no internal data found for organism ", organism))
    file_suffix <- .make_intdata_file_suffix(for.aa, domain_system)
    intdata_filename <- paste0(organism, file_suffix)
    intdata_path <- file.path(intdata_dir, intdata_filename)
    if (!file.exists(intdata_path))
        stop(wmsg("internal data file ", intdata_filename, " ",
                  "not found in ", intdata_dir))
    intdata_path
}

make_germline_db_intdata_path <- function(db_path, for.aa, domain_system)
{
    stopifnot(dir.exists(db_path))
    file_suffix <- .make_intdata_file_suffix(for.aa, domain_system)
    intdata_filename <- paste0("V", file_suffix)
    file.path(db_path, "internal_data", intdata_filename)
}

get_intdata_path <- function(organism, for.aa=FALSE,
                             domain_system=c("imgt", "kabat"),
                             which=c("live", "original"))
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    if (!isTRUEorFALSE(for.aa))
        stop(wmsg("'for.aa' must be TRUE or FALSE"))
    domain_system <- match.arg(domain_system)
    which <- match.arg(which)

    if (!valid_germline_db_name(organism))
        return(.get_igblast_intdata_path(which, organism,
                                         for.aa, domain_system))

    ## Treat 'organism' as a valid germline db name.
    db_path <- get_germline_db_path(organism)
    intdata_path <- make_germline_db_intdata_path(db_path,
                                                  for.aa, domain_system)
    if (!dir.exists(dirname(intdata_path)))
        stop(wmsg("no internal data found in germline db ", organism))
    if (!file.exists(intdata_path))
        stop(wmsg("internal data file ", basename(intdata_path), " ",
                  "not found in germline db ", organism))
    intdata_path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_intdata()
###

### IMPORTANT NOTE: The FWR/CDR positions in the returned data.frame are
### 1-based while the coding frame start positions are 0-based!
load_intdata <- function(organism, for.aa=FALSE,
                         domain_system=c("imgt", "kabat"),
                         which=c("live", "original"))
{
    domain_system <- match.arg(domain_system)
    which <- match.arg(which)
    intdata_path <- get_intdata_path(organism, for.aa=for.aa,
                                     domain_system=domain_system,
                                     which=which)
    read_ndm_data(intdata_path)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### V_genes_with_varying_fwrcdr_boundaries()
###

.get_intdata_col <- function(intdata, colname)
{
    if (!is.data.frame(intdata))
        stop(wmsg("'intdata' must be a data.frame as returned ",
                  "by load_intdata()"))
    if (!isSingleNonWhiteString(colname))
        stop(wmsg("'colname' must be a single (non-empty) string"))
    intdata_col <- intdata[[colname]]
    if (is.null(intdata_col))
        stop(wmsg("'intdata' has no \"", colname, "\" column. Make sure ",
                  "that it's a data.frame as returned by load_intdata()."))
    intdata_col
}

.extract_gene_names_as_factor <- function(intdata)
{
    allele_names <- .get_intdata_col(intdata, "allele_name")
    gene_names <- allele2gene(allele_names)
    unique_gene_names <- unique(gene_names)
    factor(gene_names, levels=unique_gene_names)
}

.check_V_segment <- function(V_segment)
{
    if (!isSingleNonWhiteString(V_segment))
        stop(wmsg("'V_segment' must be a single (non-empty) string"))
    if (!(V_segment %in% V_GENE_SEGMENTS)) {
        in1string <- paste0("\"", V_GENE_SEGMENTS, "\"", collapse=", ")
        stop(wmsg("'V_segment' must be one of ", in1string))
    }
}

.V_genes_with_varying_segment_boundaries <- function(intdata, V_segment)
{
    f <- .extract_gene_names_as_factor(intdata)
    .check_V_segment(V_segment)
    starts <- .get_intdata_col(intdata, paste0(V_segment, "_start"))
    ends <- .get_intdata_col(intdata, paste0(V_segment, "_end"))
    starts_per_gene <- unique(splitAsList(starts, f))
    ends_per_gene <- unique(splitAsList(ends, f))
    levels(f)[lengths(starts_per_gene) != 1L | lengths(ends_per_gene) != 1L]
}

V_genes_with_varying_fwrcdr_boundaries <- function(intdata, V_segment=NULL)
{
    if (!is.null(V_segment))
        return(.V_genes_with_varying_segment_boundaries(intdata, V_segment))
    found_genes <- lapply(V_GENE_SEGMENTS,
        function(V_segment)
            .V_genes_with_varying_segment_boundaries(intdata, V_segment))
    found_genes <- unique(unlist(found_genes, use.names=FALSE))
    ## Return the gene names in the same order as they show up in 'intdata'.
    unique_gene_names <- levels(.extract_gene_names_as_factor(intdata))
    unique_gene_names[unique_gene_names %in% found_genes]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_V_alleles()
### V_allele_has_stop_codon()
###

### Extracts the specified column from the 'indata' data.frame, and
### subset/reorder it to keep only the column values that correspond
### to the alleles in 'V_alleles'. Returns them in a named vector that
### is parallel to 'V_alleles' and has the allele names on it.
### The returned vector will have NAs for alleles that are not annotated
### in 'indata' or when 'indata[[colname]]' reports an NA for the allele.
.query_intdata <- function(intdata, V_alleles, colname)
{
    allele_names <- .get_intdata_col(intdata, "allele_name")
    if (!is(V_alleles, "DNAStringSet"))
        stop(wmsg("'V_alleles' must be DNAStringSet object"))
    V_names <- names(V_alleles)
    if (is.null(V_names))
        stop(wmsg("'V_alleles' must have names"))
    intdata_col <- .get_intdata_col(intdata, colname)
    setNames(intdata_col[match(V_names, allele_names)], V_names)
}

.translate_V_codons <- function(V_alleles, offsets, with.init.codon)
{
    stopifnot(is(V_alleles, "DNAStringSet"), is.integer(offsets),
              length(V_alleles) == length(offsets))
    ans <- rep.int(NA_character_, length(V_alleles))
    selection_idx <- which(!is.na(offsets))
    if (length(selection_idx) != 0L) {
        dna <- V_alleles[selection_idx]
        off <- offsets[selection_idx]
        aa <- translate_codons(dna, offset=off, with.init.codon=with.init.codon)
        ans[selection_idx] <- as.character(aa)
    }
    setNames(ans, names(V_alleles))
}

### Translates the coding frame contained in the V allele sequence.
### Only needs access to the "coding_frame_start" column in 'intdata'.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'V_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'intdata' or for which 'intdata$coding_frame_start' has
### an NA.
.translate_V_coding_frame <- function(V_alleles, intdata)
{
    offsets <- .query_intdata(intdata, V_alleles, "coding_frame_start")
    .translate_V_codons(V_alleles, offsets, with.init.codon=TRUE)
}

### Only needs access to the "<V_segment>_start" and "<V_segment>_end"
### columns of the 'intdata' data.frame.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'V_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is
### not annotated in 'intdata' or for which 'intdata$<V_segment>_start'
### or 'intdata$<V_segment>_end' has an NA.
.translate_V_segment <- function(V_alleles, intdata, V_segment)
{
    .check_V_segment(V_segment)
    start_colname <- paste0(V_segment, "_start")
    end_colname <- paste0(V_segment, "_end")
    starts <- .query_intdata(intdata, V_alleles, start_colname)  # 1-based
    ends <- .query_intdata(intdata, V_alleles, end_colname)  # 1-based
    offsets <- starts - 1L
    with.init.codon <- V_segment == "fwr1"
    ans <- .translate_V_codons(V_alleles, offsets, with.init.codon)
    ncodons <- (ends - offsets) %/% 3L
    substr(ans, 1L, ncodons)
}

translate_V_alleles <- function(V_alleles, intdata, V_segment=NULL)
{
    if (is.null(V_segment))
        return(.translate_V_coding_frame(V_alleles, intdata))
    .translate_V_segment(V_alleles, intdata, V_segment)
}

### Only needs access to the "coding_frame_start" column in 'intdata'.
### Returns a named logical vector that is parallel to 'V_alleles' and has
### the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'intdata' or for which 'intdata$coding_frame_start' has an NA.
V_allele_has_stop_codon <- function(V_alleles, intdata)
{
    V_aa <- translate_V_alleles(V_alleles, intdata)
    ans <- setNames(grepl("*", V_aa, fixed=TRUE), names(V_aa))
    ans[is.na(V_aa)] <- NA
    ans
}


### =========================================================================
### Access IgBLAST internal data
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

.extract_gene_names_as_factor <- function(intdata)
{
    allele_names <- get_intdata_col(intdata, "allele_name")
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
    starts <- get_intdata_col(intdata, paste0(V_segment, "_start"))
    ends <- get_intdata_col(intdata, paste0(V_segment, "_end"))
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


### =========================================================================
### Miscellaneous internal data utilities
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### write_intdata_to_db()
###

### Not exported!
### Note that we don't handle custom "pdm" data at the moment. So the
### function does not need a 'for.aa' argument. Also the name of the
### function reflects the fact that it only handles "ndm" data.
write_intdata_to_db <- function(ndm_data, db_path,
                                domain_system=c("imgt", "kabat"),
                                check.and.reorder=FALSE)
{
    stopifnot(is.data.frame(ndm_data), isTRUEorFALSE(check.and.reorder))
    domain_system <- match.arg(domain_system)
    intdata_path <- make_germline_db_intdata_path(db_path, FALSE, domain_system)
    intdata_dir <- dirname(intdata_path)
    stopifnot(!dir.exists(intdata_dir))

    ## Even though write_ndm_data() will call check_ndm_data_col2class()
    ## internally, we prefer to fail **before** creating the 'intdata_dir'
    ## folder.
    check_ndm_data_col2class(ndm_data)
    if (check.and.reorder) {
        db_V_fasta_file <- get_db_fasta_file(db_path, "V")
        db_V_allele_names <- names(fasta.seqlengths(db_V_fasta_file))
        ndm_data <- check_and_reorder_igdata_rows(ndm_data, db_V_allele_names)
    }
    stopifnot(dir.create(intdata_dir))
    write_ndm_data(ndm_data, intdata_path)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show_intdata_disagreements()
###

show_intdata_disagreements <- function(db_name)
{
    check_germline_db_name(db_name)
    db_intdata <- load_intdata(db_name)
    igblast_organism <- infer_igblast_organism_from_db_name(db_name)
    if (is.na(igblast_organism))
        stop(wmsg("The specified germline db does not seem to be ",
                  "for an IgBLAST organism. Note that you can use ",
                  "list_igblast_organisms() to get the list of ",
                  "IgBLAST organisms. See '?list_igblast_organisms' ",
                  "for more information."))
    igblast_intdata <- load_intdata(igblast_organism)
    diff <- df_diff(db_intdata, igblast_intdata, "allele_name",
                    "igblastr-generated", "IgBLAST-provided")
    what <- c("the igblastr-generated \"internal data\" included in this ",
              "germline db and the \"internal data\" provided by IgBLAST ",
              "for ", igblast_organism)
    if (length(diff) == 0L) {
        msg <- c("No disagreements between ", what, ".")
        cat(wmsg2(msg, margin=0L), "\n", sep="")
    } else {
        msg <- c("Disagreements between ", what, ":")
        cat(wmsg2(msg, margin=0L), "\n", sep="")
        cat(diff, sep="")
    }
}


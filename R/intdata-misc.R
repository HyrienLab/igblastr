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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### count_conserved_V_codons()
### summarize_anchor_V_residues()
###

.get_amino_acids_at <- function(V_alleles, codon_ends)
{
    stopifnot(is(V_alleles, "DNAStringSet"), is.integer(codon_ends),
              length(V_alleles) == length(codon_ends))
    codon_starts <- codon_ends - 2L
    idx <- which(codon_starts >= 1L & codon_ends <= width(V_alleles))
    codons <- subseq(V_alleles[idx], codon_starts[idx], codon_ends[idx])
    translate(codons, no.init.codon=TRUE, if.fuzzy.codon="solve")
}

### Not exported!
### If the 2nd argument ('intdata') is missing then the 1st argument
### ('V_alleles') must be the name of an existing germline db.
### IMGT codons 23 (fwr1), 41 (fwr2), and 104 (fwr3) are expected to
### translate to C, W, and C, respectively (conserved codons).
### Returns a data.frame with 1 row per codon of interest and the
### following columns: C_count, W_count, other, percent_C, percent_W.
count_conserved_V_codons <- function(V_alleles, intdata)
{
    if (missing(intdata)) {
        check_germline_db_name(V_alleles)
        intdata <- load_intdata(V_alleles)
        V_alleles <- load_germline_sequences(V_alleles, region_types="V")
    }
    fwr1_ends   <- query_intdata(intdata, V_alleles, "fwr1_end", no.NAs=TRUE)
    fwr2_starts <- query_intdata(intdata, V_alleles, "fwr2_start", no.NAs=TRUE)
    fwr3_ends   <- query_intdata(intdata, V_alleles, "fwr3_end", no.NAs=TRUE)
    codons <- list(
        ## codon 23 + a few flanking codons:
        codon22 =.get_amino_acids_at(V_alleles, fwr1_ends - 12L),
        codon23 =.get_amino_acids_at(V_alleles, fwr1_ends - 9L),
        codon24 =.get_amino_acids_at(V_alleles, fwr1_ends - 6L),
        codon25 =.get_amino_acids_at(V_alleles, fwr1_ends - 3L),
        codon26 =.get_amino_acids_at(V_alleles, fwr1_ends),
        ## codon 41 + a few flanking codons:
        codon39 =.get_amino_acids_at(V_alleles, fwr2_starts + 2L),
        codon40 =.get_amino_acids_at(V_alleles, fwr2_starts + 5L),
        codon41 =.get_amino_acids_at(V_alleles, fwr2_starts + 8L),
        codon42 =.get_amino_acids_at(V_alleles, fwr2_starts + 11L),
        codon43 =.get_amino_acids_at(V_alleles, fwr2_starts + 14L),
        ## codon 104 + a few flanking codons:
        codon103=.get_amino_acids_at(V_alleles, fwr3_ends - 3L),
        codon104=.get_amino_acids_at(V_alleles, fwr3_ends),
        codon105=.get_amino_acids_at(V_alleles, fwr3_ends + 3L)
    )
    C_count <- vapply(codons, function(codon) sum(codon == "C"), integer(1))
    W_count <- vapply(codons, function(codon) sum(codon == "W"), integer(1))
    other <- vapply(codons, function(codon) sum(codon != "C" & codon != "W"),
                    integer(1))
    percent_C <- round(100 * C_count / lengths(codons), digits=2L)
    percent_W <- round(100 * W_count / lengths(codons), digits=2L)
    ans <- data.frame(C_count=C_count, W_count=W_count, other=other,
                      percent_C=percent_C, percent_W=percent_W)
    rownames(ans) <- names(codons)
    ans
}

### Not exported!
### If the 2nd argument ('intdata') is missing then the 1st argument
### ('V_alleles') must be the name of an existing germline db.
### Returns C23/W41/C104 percentages in a named numeric vector.
summarize_anchor_V_residues <- function(V_alleles, intdata)
{
    conserved_codons <- c(codon23="C", codon41="W", codon104="C")
    Mindex <- cbind(names(conserved_codons),
                    paste0("percent_", conserved_codons))
    ans <- count_conserved_V_codons(V_alleles, intdata)[Mindex]
    names(ans) <- paste0(conserved_codons,
                         sub("^codon", "", names(conserved_codons)))
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### summarize_anchor_V_residues_for_IMGT_organisms()
###

### Not exported!
### Percent cysteines (C) at codon23 & codon104 and percent tryptophans (W)
### at codon41. Results obtained with igblastr 1.3.20 for IMGT 202614-2:
### ---------------------------- tcr.db=FALSE --------------------------------
###                                                          C23    W41   C104
### IMGT-202614-2.Bos_taurus.IGH+IGK+IGL                   98.88  95.51  95.00
### IMGT-202614-2.Canis_lupus_familiaris.IGH+IGK+IGL       97.89  47.68  47.68
### IMGT-202614-2.Equus_caballus.IGH+IGK                   91.01  88.76  34.83
### IMGT-202614-2.Gallus_gallus.IGH+IGL                    73.98  67.91  71.09
### IMGT-202614-2.Gorilla_gorilla_gorilla.IGH+IGK+IGL      98.26  89.58  93.40
### IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL                 98.47  94.52  97.41
### IMGT-202614-2.Lemur_catta.IGH+IGK+IGL                  96.92  99.56 100.00
### IMGT-202614-2.Macaca_fascicularis.IGH                  93.18   0.00   0.00
### IMGT-202614-2.Macaca_mulatta.IGH+IGK+IGL               95.62  93.87  96.05
### IMGT-202614-2.Mus_musculus.IGH+IGK+IGL                 98.49  94.68  94.51
### IMGT-202614-2.Mustela_putorius_furo.IGH+IGK+IGL        65.54  31.64  32.20
### IMGT-202614-2.Neogale_vison.IGH                         0.00   0.00   0.00
### IMGT-202614-2.Oncorhynchus_mykiss.IGH                  98.73  92.36  94.90
### IMGT-202614-2.Ornithorhynchus_anatinus.IGH             97.78   0.00   0.00
### IMGT-202614-2.Oryctolagus_cuniculus.IGH+IGK+IGL        98.64 100.00 100.00
### IMGT-202614-2.Pongo_pygmaeus.IGH+IGK+IGL               49.66  55.48  48.11
### IMGT-202614-2.Rattus_norvegicus.IGH+IGK+IGL            95.76  40.20  42.21
### IMGT-202614-2.Salmo_salar.IGH                           0.00   0.00   0.67
### IMGT-202614-2.Sus_scrofa.IGH+IGK+IGL                   96.92  96.92  95.38
### IMGT-202614-2.Vicugna_pacos.IGH                        98.81  97.62 100.00
### ---------------------------- tcr.db=TRUE ---------------------------------
###                                                          C23    W41   C104
### IMGT-202614-2.Bos_taurus.TRA+TRB+TRG+TRD               66.58  48.29  48.56
### IMGT-202614-2.Camelus_dromedarius.TRA+TRB+TRG+TRD     100.00 100.00 100.00
### IMGT-202614-2.Canis_lupus_familiaris.TRA+TRB+TRG+TRD   98.77  98.77  14.81
### IMGT-202614-2.Danio_rerio.TRA+TRD                       0.71   0.00   0.00
### IMGT-202614-2.Felis_catus.TRA+TRB+TRG+TRD              98.85  98.85  44.83
### IMGT-202614-2.Gorilla_gorilla_gorilla.TRA+TRB+TRG+TRD  53.72  54.79  45.21
### IMGT-202614-2.Heterocephalus_glaber.TRA+TRB+TRG+TRD    95.65  95.65  90.43
### IMGT-202614-2.Homo_sapiens.TRA+TRB+TRG+TRD             96.88  98.87  97.69
### IMGT-202614-2.Macaca_fascicularis.TRB                  95.45  98.48  98.46
### IMGT-202614-2.Macaca_mulatta.TRA+TRB+TRG+TRD           34.26  34.72  62.79
### IMGT-202614-2.Mus_musculus.TRA+TRB+TRG+TRD             99.47  97.16  98.70
### IMGT-202614-2.Mustela_putorius_furo.TRA+TRB+TRG+TRD    96.64  94.12  36.97
### IMGT-202614-2.Oryctolagus_cuniculus.TRA+TRB+TRG+TRD    99.32  99.32 100.00
### IMGT-202614-2.Ovis_aries.TRA+TRB+TRD                   99.47  72.80  18.72
### IMGT-202614-2.Pan_troglodytes.TRA+TRG+TRD              97.50  98.75  81.25
### IMGT-202614-2.Pongo_abelii.TRA+TRB+TRG+TRD             95.83  96.67  85.00
### IMGT-202614-2.Pongo_pygmaeus.TRB+TRG                   95.73  96.58  84.62
### IMGT-202614-2.Sus_scrofa.TRB+TRG                       92.68  97.56   0.00
summarize_anchor_V_residues_for_IMGT_organisms <-
    function(release, tcr.db=FALSE)
{
    imgt_organisms <- list_IMGT_organisms(release)
    all_percents <- lapply(imgt_organisms,
        function(organism) {
            message(organism)
            db_name <- try(suppressWarnings(suppressMessages(
                install_IMGT_germline_db(release, organism,
                                         tcr.db=tcr.db, overwrite=TRUE)
            )), silent=TRUE)
            if (inherits(db_name, "try-error"))
                return(NULL)
            percents <- summarize_anchor_V_residues(db_name)
            matrix(percents, nrow=1L, dimnames=list(db_name, names(percents)))
        })
    do.call(rbind, all_percents)
}


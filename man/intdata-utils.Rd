\name{intdata-utils}

\alias{intdata-utils}
\alias{intdata_utils}
\alias{intdata}

\alias{get_intdata_path}
\alias{load_intdata}
\alias{V_genes_with_varying_fwrcdr_boundaries}
\alias{translate_V_alleles}
\alias{V_allele_has_stop_codon}
\alias{annotate_heavy_V_alleles}
\alias{annotate_light_V_alleles}

\title{Access and manipulate IgBLAST internal data}

\description{
  IgBLAST \emph{internal data} is expected to annotate all the known
  germline V gene alleles for a given organism. It is provided by NCBI
  and is typically included in a standard IgBLAST installation.

  The \pkg{igblastr} package provides a small set of utilities to access
  and manipulate IgBLAST \emph{internal data}.
}

\usage{
## Access internal data:
get_intdata_path(organism, for.aa=FALSE, domain_system=c("imgt", "kabat"),
                 which=c("live", "original"))
load_intdata(organism, for.aa=FALSE, domain_system=c("imgt", "kabat"),
             which=c("live", "original"))

## Manipulate internal data:
V_genes_with_varying_fwrcdr_boundaries(intdata, V_segment=NULL)
translate_V_alleles(V_alleles, intdata, V_segment=NULL)
V_allele_has_stop_codon(V_alleles, intdata)
}

\arguments{
  \item{organism}{
    A single string containing the name of an organism as
    returned by \code{\link{list_igblast_organisms}()}.

    Alternatively, this can be the name of a cached germline db.
    Note that this works only for germline dbs that include their
    own \emph{internal data}. You can use:
    \preformatted{    list_germline_dbs(with.intdata.only=TRUE)}
    to list them.
    See \code{?\link{list_germline_dbs}} for more information.
  }
  \item{for.aa}{
    By default, the data.frame returned by \code{load_intdata()} contains
    FWR/CDR boundaries reported with respect to the nucleotide sequences
    of the germline V alleles.
    Setting \code{for.aa} to \code{TRUE} will return a data.frame where
    they are reported with respect to the amino acid sequences of the
    germline V alleles.
  }
  \item{domain_system}{
    Domain system to be used for segment annotation. Must be \code{"imgt"}
    (the default) or \code{"kabat"}.
  }
  \item{which}{
    By default, \code{get_intdata_path()} and \code{load_intdata()}
    access the "live IgBLAST data", that is, the IgBLAST data that
    the user has possibly updated with \code{update_live_igdata()}.
    Depending on whether updates were applied or not, the "live IgBLAST data"
    might differ from the original IgBLAST data.

    Set \code{which} to \code{"original"} if you want to access
    the original IgBLAST data instead.

    See \code{?\link{update_live_igdata}} for more information about
    "live" and "original" IgBLAST data.
  }
  \item{intdata}{
    A data.frame as returned by \code{load_intdata()} for
    \code{V_genes_with_varying_fwrcdr_boundaries()}, or by
    \code{load_intdata(..., for.aa=FALSE)} for \code{translate_V_alleles()}
    and \code{V_allele_has_stop_codon()}.
  }
  \item{V_segment}{
    The name of a V gene segment. This can be set to \code{"fwr1"},
    \code{"cdr1"}, \code{"fwr2"}, \code{"cdr2"}, or \code{"fwr3"}.

    By default (i.e. when \code{V_segment} is omitted or set to \code{NULL}),
    \code{V_genes_with_varying_fwrcdr_boundaries()} will identify V genes for
    which any segment has varying boundaries across alleles. Otherwise, it
    will identify V genes for which the specified segment has varying
    boundaries.

    By default \code{translate_V_alleles()} will translate the entire coding
    frame in each allele. Otherwise, it will translate the specified segment
    only.
  }
  \item{V_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline V
    gene allele sequences.
  }
}

\details{
  IgBLAST \emph{internal data} is typically included in a standard IgBLAST
  installation. It's located in the \code{internal_data/} directory which
  is itself a subdirectory of IgBLAST \emph{root directory}.
}

\value{
  \code{get_intdata_path()} returns a single string containing
  the path to the \emph{internal data} included in the IgBLAST
  installation used by \pkg{igblastr}, for the specified organism.

  \code{load_intdata()} returns the \emph{internal data} in a data.frame
  with 1 row per germline V gene allele sequence and the same columns
  as the data.frame returned by \code{read_ndm_data()}.
  See \code{?\link{read_ndm_data}} for more information.

  \code{V_genes_with_varying_fwrcdr_boundaries()} returns a character
  vector containing the names of the germline V genes for which the
  FWR/CDR boundaries are not the same across all alleles.

  \code{translate_V_alleles()} returns a named character vector with 1
  amino acid sequence per supplied allele. The vector contains an \code{NA}
  for any allele that is not annotated in \code{intdata} or for which
  the required information is \code{NA}. The names on it are
  the names of the supplied alleles.

  \code{V_allele_has_stop_codon()} returns a named logical vector with 1
  value per supplied allele. The vector contains an \code{NA} for any
  allele that is not annotated in \code{intdata} or for which
  \code{intdata$coding_frame_start} has an \code{NA}. The names on it are
  the names of the supplied alleles.
}

\seealso{
  \itemize{
    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.

    \item \link{auxdata_utils} to access, manipulate, and generate IgBLAST
          \emph{auxiliary data}.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \code{\link{list_igblast_organisms}} to list the organisms
          supported by IgBLAST.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link{translate_codons}} function on which
          \code{translate_V_alleles()} is based.

    \item \code{\link{allele2gene}} to go from germline gene allele names
          to germline gene names.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## list_igblast_organisms() and get_intdata_path()
## ---------------------------------------------------------------------

list_igblast_organisms()

get_intdata_path("rabbit")
rabbit_intdata <- load_intdata("rabbit")
head(rabbit_intdata)

rabbit_intdata2 <- load_intdata("rabbit", for.aa=TRUE)
head(rabbit_intdata2)

## The values in the "end" cols in 'rabbit_intdata' are exactly 3 times
## those in the "end" cols in 'rabbit_intdata2':
end_colnames <- grep("_end$", colnames(rabbit_intdata), value=TRUE)
stopifnot(identical(rabbit_intdata [ , end_colnames],
                    rabbit_intdata2[ , end_colnames] * 3L))

## Get the internal data included germline db
## _OGRDB.human.IGH+IGK+IGL.202410 (this data
## is provided by AIRR-community/OGRDB):
db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
human_intdata <- load_intdata(db_name)
head(human_intdata)

## ---------------------------------------------------------------------
## V_genes_with_varying_fwrcdr_boundaries()
## ---------------------------------------------------------------------

## Note that the alleles of a given germline V gene don't necessarily
## share the same FWR/CDR boundaries. You can use utility function
## V_genes_with_varying_fwrcdr_boundaries() to identify them:
human_intdata0 <- load_intdata("human")
var_genes <- V_genes_with_varying_fwrcdr_boundaries(human_intdata0)
var_genes
subset(human_intdata0, allele2gene(allele_name) == "IGHV4-31")

## Human germline V genes for which the CDR1 boundaries are not the same
## across all alleles:
var_genes <- V_genes_with_varying_fwrcdr_boundaries(human_intdata0,
                                                    V_segment="cdr1")
var_genes
subset(human_intdata0, allele2gene(allele_name) \%in\% var_genes)

## ---------------------------------------------------------------------
## translate_V_alleles() and V_allele_has_stop_codon()
## ---------------------------------------------------------------------

V_alleles <- load_germline_db(db_name, region_types="V")
V_alleles  # DNAStringSet object

V_aa <- translate_V_alleles(V_alleles, human_intdata)
head(V_aa)

fwr2 <- translate_V_alleles(V_alleles, human_intdata, V_segment="fwr2")
head(fwr2)

## Surprisingly, 13 V alleles in _OGRDB.human.IGH+IGK+IGL.202410 contain
## the stop codon:
has_stop_codon <- grepl("*", V_aa, fixed=TRUE)
table(has_stop_codon)
V_aa[has_stop_codon]
V_alleles[has_stop_codon]
}

\keyword{utilities}

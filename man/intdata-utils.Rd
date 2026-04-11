\name{intdata-utils}

\alias{intdata-utils}
\alias{intdata_utils}
\alias{intdata}

\alias{get_intdata_path}
\alias{load_intdata}
\alias{V_genes_with_varying_fwrcdr_boundaries}

\title{Access IgBLAST internal data}

\description{
  IgBLAST \emph{internal data} is expected to annotate all the known
  germline V gene alleles for a given organism. It is provided by NCBI
  and is typically included in a standard IgBLAST installation.

  The \emph{internal data} informs about the FWR/CDR boundaries on
  the V allele sequences, also known as the V gene delineations.
  This information is needed and used internally by IgBLAST.

  \code{get_intdata_path()} and \code{load_intdata()} can be used
  to access the \emph{internal data} included in IgBLAST or in
  one of the \emph{cached germline dbs} managed by \pkg{igblastr}.

  \code{V_genes_with_varying_fwrcdr_boundaries()} can be used to identify
  germline V genes for which the FWR/CDR boundaries are not the same across
  all alleles.
}

\usage{
## Access internal data:
get_intdata_path(organism, for.aa=FALSE, domain_system=c("imgt", "kabat"),
                 which=c("live", "original"))
load_intdata(organism, for.aa=FALSE, domain_system=c("imgt", "kabat"),
             which=c("live", "original"))

## Identify V genes with varying FWR/CDR boundaries across alleles:
V_genes_with_varying_fwrcdr_boundaries(intdata, V_segment=NULL)
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
    A data.frame as returned by \code{load_intdata()}.
  }
  \item{V_segment}{
    The name of a V gene segment. This can be set to \code{"fwr1"},
    \code{"cdr1"}, \code{"fwr2"}, \code{"cdr2"}, or \code{"fwr3"}.

    By default (i.e. when \code{V_segment} is omitted or set to \code{NULL}),
    \code{V_genes_with_varying_fwrcdr_boundaries()} will identify V genes for
    which any segment has varying boundaries across alleles. Otherwise, it
    will identify V genes for which the specified segment has varying
    boundaries.
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
}

\seealso{
  \itemize{
    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.

    \item \link{auxdata_utils} to access IgBLAST \emph{auxiliary data}.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \code{\link{list_igblast_organisms}} to list the organisms
          supported by IgBLAST.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

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
## Access IgBLAST internal data for a given organism
## ---------------------------------------------------------------------

## IgBLAST only includes internal data for the following organisms:
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

## ---------------------------------------------------------------------
## Access the internal data included in a germline db
## ---------------------------------------------------------------------

## List germline dbs with internal data:
list_germline_dbs(with.intdata.only=TRUE)

## Access the internal data included in germline db
## _OGRDB.human.IGH+IGK+IGL.202410 (this data was inferred from
## the gaps in the germline V gene allele sequences provided by
## AIRR-community/OGRDB):
db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"

get_intdata_path(db_name)

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
var_genes  # 4 genes

## Display rows associated with gene IGHV4-31:
subset(human_intdata0, allele2gene(allele_name) == "IGHV4-31")

## Human germline V genes for which the CDR1 boundaries are not the
## same across all alleles:
var_genes <- V_genes_with_varying_fwrcdr_boundaries(human_intdata0,
                                                    V_segment="cdr1")
var_genes  # 2 genes

## Display rows associated with the genes in 'var_genes':
subset(human_intdata0, allele2gene(allele_name) \%in\% var_genes)
}

\keyword{utilities}

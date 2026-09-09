\name{print_gapped_V_alleles}

\alias{print_gapped_V_alleles}

\title{Fancy display of V alleles sequences}

\description{
  Use \code{print_gapped_V_alleles()} to display a set of gapped
  V allele sequences where the sequences are colored and split by
  FWR/CDR regions.
}

\usage{
print_gapped_V_alleles(gapped_V_alleles,
                       fwrcdr_widths=IMGT_FWRCDR_WIDTHS,
                       translate=FALSE,
                       region_types=NULL,
                       igblast_organism=NA,
                       filler=".", sep=" ")
}

\arguments{
  \item{gapped_V_alleles}{
    COMING SOON...
  }
  \item{fwrcdr_widths}{
    COMING SOON...
  }
  \item{translate}{
    COMING SOON...
  }
  \item{region_types}{
    COMING SOON...
  }
  \item{igblast_organism}{
    COMING SOON...
  }
  \item{filler}{
    COMING SOON...
  }
  \item{sep}{
    COMING SOON...
  }
}

\value{
  COMING SOON...
}

\seealso{
  \itemize{
    \item \link{V_alleles_inspect} for basic inspection of V allele sequences.

    \item \code{\link{load_intdata}} to access IgBLAST \emph{internal data}.

    \item \code{\link{print_J_alleles}} to display a set of J allele
          sequences where the sequences are justified w.r.t. their
          CDR3/FWR4 junction.

    \item \link{J_alleles_inspect} for basic inspection of J allele sequences.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \link[Biostrings]{DNAStringSet} and \link[Biostrings]{AAStringSet}
          objects in the \pkg{Biostrings} package.

    \item The \code{\link{translate_codons}} function which is used internally
          by \code{print_gapped_V_alleles()} when \code{translate=TRUE}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## COMING SOON...
}

\keyword{utilities}

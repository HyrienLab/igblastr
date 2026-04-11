\name{allele2gene}

\alias{allele2gene}

\title{Go from germline gene allele names to germline gene names}

\description{
  A simple convenience function to remove the allele suffix from a vector
  of germline gene allele names.
}

\usage{
allele2gene(allele_names)
}

\arguments{
  \item{allele_names}{
    A character vector of germline gene allele names.
  }
}

\details{
  Germline gene allele names use specific nomenclature, often including
  a gene name followed by an asterisk and allele number, like in
  IGHD3-16*03.

  The \code{allele2gene()} function simply removes the asterisk
  and anything that follows it to keep the gene name only. Note that
  the function is \emph{vectorized}, that is, its input can be a
  \emph{vector} of germline gene allele names, in which case the
  corresponding germline gene names are returned in a vector of
  the same length as the input vector.
  The names on the input vector are propagated, and so are the
  \code{NA}'s in it.
}

\value{
  A character vector \emph{parallel} to the input vector.
}

\seealso{
  \itemize{
    \item \code{\link{load_germline_db}} to load the nucleotide sequences
          of the gene regions stored in a cached germline db.

    \item \code{\link{load_c_region_db}} to load the nucleotide sequences of
          the gene regions stored in a cached C-region db.

    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.
  }
}

\examples{
allele2gene(c("IGHV1-2*04", "IGHV1-2*06", "IGHV5-51*01", "IGHD3-16*03"))

J_alleles <- load_germline_db("_OGRDB.human.IGH+IGK+IGL.202410",
                              region_types="J")
names(J_alleles)
allele2gene(names(J_alleles))

C_alleles <- load_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")
names(C_alleles)
allele2gene(names(C_alleles))
}

\keyword{utilities}

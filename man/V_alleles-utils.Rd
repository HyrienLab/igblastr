\name{V_alleles-utils}

\alias{V_alleles-utils}
\alias{V_alleles_utils}

\alias{translate_V_alleles}
\alias{V_allele_has_stop_codon}

\title{Basic manipulation of V allele sequences}

\description{
  A small set of utilities for basic examination of V allele sequences.

  Note that all of them require access to the \emph{internal data}
  associated with the V alleles to examine. See \code{intdata} argument
  below for how to obtain this data.
}

\usage{
translate_V_alleles(V_alleles, intdata, V_segment=NULL)
V_allele_has_stop_codon(V_alleles, intdata)
}

\arguments{
  \item{V_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline V
    gene allele sequences.
  }
  \item{intdata}{
    A data.frame as returned by \code{\link{load_intdata}(..., for.aa=FALSE)}.
  }
  \item{V_segment}{
    The name of a V gene segment. This can be set to \code{"fwr1"},
    \code{"cdr1"}, \code{"fwr2"}, \code{"cdr2"}, or \code{"fwr3"}.

    By default \code{translate_V_alleles()} will translate the entire coding
    frame in each allele. Otherwise, it will translate the specified segment
    only.
  }
}

\value{
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
    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

    \item \link{auxdata_utils} to access IgBLAST \emph{auxiliary data}.

    \item \link{J_alleles_utils} for basic manipulation of J allele sequences.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link{translate_codons}} function upon which
          \code{translate_V_alleles()} is based.

    \item \code{\link{allele2gene}} to go from germline gene allele names
          to germline gene names.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
## We're going to examine the V allele sequences stored in
## germline db _OGRDB.human.IGH+IGK+IGL.202410.

## Load the V allele sequences:
db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
V_alleles <- load_germline_db(db_name, region_types="V")
V_alleles  # DNAStringSet object

## Load the corresponding internal data:
intdata <- load_intdata(db_name)

## Translate the V allele sequences:
V_aa <- translate_V_alleles(V_alleles, intdata)
head(V_aa)

## Translate the FWR2 portion of the V allele sequences:
fwr2_aa <- translate_V_alleles(V_alleles, intdata, V_segment="fwr2")
head(fwr2_aa)

## Surprisingly, 13 V alleles in _OGRDB.human.IGH+IGK+IGL.202410
## contain the stop codon:
has_stop_codon <- grepl("*", V_aa, fixed=TRUE)
table(has_stop_codon)

## Display the V allele sequences with stop codon:
V_aa[has_stop_codon]       # protein
V_alleles[has_stop_codon]  # nucleotides
}

\keyword{utilities}

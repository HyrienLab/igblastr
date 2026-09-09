\name{J_alleles-inspect}

\alias{J_alleles-inspect}
\alias{J_alleles_inspect}

\alias{translate_J_alleles}
\alias{J_allele_has_stop_codon}
\alias{translate_fwr4}

\title{Basic inspection of J allele sequences}

\description{
  A small set of utilities for basic inspection of J allele sequences.

  Note that all of the functions documented in this man page require
  access to the \emph{auxiliary data} associated with the J alleles to
  inspect. See \code{auxdata} argument below for how to obtain this data.
}

\usage{
translate_J_alleles(J_alleles, auxdata)
J_allele_has_stop_codon(J_alleles, auxdata)
translate_fwr4(J_alleles, auxdata, max.codons=NA)
}

\arguments{
  \item{J_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline J
    gene allele sequences.
  }
  \item{auxdata}{
    A data.frame as returned by \code{\link{load_auxdata}()},
    \code{\link{compute_auxdata}()}, or
    \code{\link{compute_germline_db_auxdata}()}, that contains
    annotations for the alleles in \code{J_alleles}.
  }
  \item{max.codons}{
    The maximum number of FWR4 codons to translate. By default (i.e. when
    \code{max.codons} is \code{NA}) all the FWR4 codons are translated.
  }
}

\value{
  \code{translate_J_alleles()} returns a named character vector with 1
  amino acid sequence per supplied allele. The vector contains an \code{NA}
  for any allele that is not annotated in \code{auxdata} or for which
  \code{auxdata$coding_frame_start} has an \code{NA}. The names on it are
  the names of the supplied alleles.

  \code{J_allele_has_stop_codon()} returns a named logical vector with 1
  value per supplied allele. The vector contains an \code{NA} for any
  allele that is not annotated in \code{auxdata} or for which
  \code{auxdata$coding_frame_start} has an \code{NA}. The names on it are
  the names of the supplied alleles.

  \code{translate_fwr4()} returns a named character vector with 1 amino
  acid sequence per supplied allele. The vector contains an \code{NA}
  for any allele that is not annotated in \code{auxdata} or for which
  \code{auxdata$cdr3_end} has an \code{NA}.
}

\seealso{
  \itemize{
    \item \code{\link{print_J_alleles}} to display a set of J allele
          sequences where the sequences are justified w.r.t. their
          CDR3/FWR4 junction.

    \item \code{\link{load_auxdata}} to access IgBLAST \emph{auxiliary data}.

    \item \link{V_alleles_inspect} for basic inspection of V allele sequences.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \link[Biostrings]{DNAStringSet} and \link[Biostrings]{AAStringSet}
          objects in the \pkg{Biostrings} package.

    \item The \code{\link{translate_codons}} function which is used internally
          by \code{print_J_alleles()}, \code{translate_J_alleles()}, and
          \code{translate_fwr4()}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## 1. translate_J_alleles()
## ---------------------------------------------------------------------

db_name <- install_IMGT_germline_db("202631-1", "Oryctolagus cuniculus",
                                    overwrite=TRUE)
J_alleles <- load_germline_sequences(db_name, region_types="J")
auxdata <- load_auxdata(db_name)
head(auxdata)

## The 'coding_frame_start' column in 'auxdata' contains integer
## values that are >= 0 and <= 2. They indicate how many nucleotides
## precede the first codon on each allele sequence. In other words,
## this is the number of nucleotides that we need to trim on the 5'
## end of the germline J allele sequence before we start translating
## it. translate_J_alleles() uses this information to translate the
## DNA sequences in 'J_alleles':
J_aa <- translate_J_alleles(J_alleles, auxdata)
J_aa

## No sequence in 'J_aa' should contain the letter "*" which is used
## by translate_J_alleles() to represent a stop codon. However, one
## J allele in IMGT-202631-1.Oryctolagus_cuniculus.IGH+IGK+IGL seems
## to disobey:
has_stop_codon <- grepl("*", J_aa, fixed=TRUE)
auxdata[has_stop_codon, ]  # coding_frame_start = 0 for IGKJ1-2*04
J_alleles[has_stop_codon]  # first codon (TGA) is a stop codon
J_aa[has_stop_codon]       # indeed!

## ---------------------------------------------------------------------
## 2. ABOUT THE "WGXG" AND "FGXG" MOTIFS
## ---------------------------------------------------------------------

## The FWR4 region is expected to start with the following amino acid
## motifs (X represents any amino acid):
##   - "WGXG" on the heavy chain
##   - "FGXG" on the light chain

## Let's use translate_fwr4() to extract and translate the first 4
## codons of the FWR4 region:
fwr4_head <- translate_fwr4(J_alleles, auxdata, max.codons=4)

## We expect to see the "WGXG" and "FGXG" motifs here, and most of the
## time we do:
has_motif <- grepl("[FW]G.G", fwr4_head)
table(has_motif)

## However, there are a few exceptions:
fwr4_head[!has_motif]
}

\keyword{utilities}

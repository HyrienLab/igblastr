\name{J_alleles-utils}

\alias{J_alleles-utils}
\alias{J_alleles_utils}

\alias{translate_J_alleles}
\alias{J_allele_has_stop_codon}
\alias{translate_fwr4}

\title{Basic manipulation of J allele sequences}

\description{
  A small set of utilities for basic examination of J allele sequences.

  Note that all of them require access to the \emph{auxiliary data}
  associated with the J alleles to examine. See \code{auxdata} argument
  below for how to obtain this data.
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
    A data.frame as returned by \code{\link{compute_auxdata}()},
    \code{\link{compute_germline_db_auxdata}()}, or
    \code{\link{load_auxdata}()}.
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
    \item \link{auxdata_utils} to access IgBLAST \emph{auxiliary data}.

    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

    \item \link{V_alleles_utils} for basic manipulation of V allele sequences.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link{translate_codons}} function upon which
          \code{translate_J_alleles()} and \code{translate_fwr4()} are based.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## 1. A close look at the J allele sequences for rabbit
## ---------------------------------------------------------------------

## IgBLAST only includes auxiliary data for the following organisms:
list_igblast_organisms()

get_auxdata_path("rabbit")
rabbit_auxdata <- load_auxdata("rabbit")

## It turns out that IgBLAST auxiliary data for rabbit matches exactly
## the set of rabbit germline J alleles available at IMGT:
db_name <- install_IMGT_germline_db("202614-2", "Oryctolagus cuniculus",
                                    without.auxdata=TRUE,
                                    overwrite=TRUE)
J_alleles <- load_germline_db(db_name, region_types="J")
J_alleles  # DNAStringSet object
stopifnot(setequal(names(J_alleles), rabbit_auxdata$allele_name))

## Note that this might change with future IMGT releases.

## Let's put the allele sequences in 'J_alleles' in the same order as
## in 'rabbit_auxdata':
J_alleles <- J_alleles[rabbit_auxdata$allele_name]
stopifnot(identical(names(J_alleles), rabbit_auxdata$allele_name))

## The 'coding_frame_start' column in 'rabbit_auxdata' contains integer
## values that are >= 0 and <= 2. They indicate how many nucleotides
## precede the first codon on each allele sequence. In other words,
## this is the number of nucleotides that we need to trim on the 5'
## end of the germline J allele sequence before we start translating
## it. translate_J_alleles() uses this information to translate the
## DNA sequences in 'J_alleles':
J_aa <- translate_J_alleles(J_alleles, rabbit_auxdata)
J_aa

## No sequence in 'J_aa' should contain the letter "*" which is used
## by translate_J_alleles() to represent a stop codon. However, one
## J allele in IMGT-202614-2.Oryctolagus_cuniculus.IGH+IGK+IGL seems
## to disobey:
has_stop_codon <- grepl("*", J_aa, fixed=TRUE)
rabbit_auxdata[has_stop_codon, ]  # coding_frame_start = 0 for IGKJ1-2*04
J_alleles[has_stop_codon]         # first codon (TGA) is a stop codon
J_aa[has_stop_codon]              # indeed!

## ---------------------------------------------------------------------
## 2. About the "WGXG" and "FGXG" motifs
## ---------------------------------------------------------------------

## The FWR4 region is expected to start with the following amino acid
## motifs (X represents any amino acid):
##   - "WGXG" on the heavy chain
##   - "FGXG" on the light chain

## Let's use translate_fwr4() to extract and translate the first 4
## codons of the FWR4 region:
fwr4_head <- translate_fwr4(J_alleles, rabbit_auxdata, max.codons=4)

## We expect to see the "WGXG" and "FGXG" motifs here, and most of the
## time we do:
has_motif <- grepl("[FW]G.G", fwr4_head)
table(has_motif)

## However, there are a few exceptions:
fwr4_head[!has_motif]

## ---------------------------------------------------------------------
## 3. Compute auxiliary data for a set of J allele sequences
## ---------------------------------------------------------------------

## See '?compute_auxdata' for the details of what compute_auxdata() does.
computed_auxdata <- compute_auxdata(J_alleles)
head(computed_auxdata)

## Alleles for which the CDR3 end could not found:
cdr3_end_not_found <- is.na(computed_auxdata$cdr3_end)
stopifnot(identical(cdr3_end_not_found, !has_motif))
J_alleles[cdr3_end_not_found]
fwr4_head[cdr3_end_not_found]  # déjà vu

## 'computed_auxdata' is in agreement with 'rabbit_auxdata', except for
## the 8 alleles for which the CDR3 end could not be found:
keep_idx <- which(!cdr3_end_not_found)
stopifnot(identical(computed_auxdata[keep_idx, ],
                    rabbit_auxdata[keep_idx, ]))
}

\keyword{utilities}

\name{translate_codons}

\alias{translate_codons}
\alias{extract_codons}
\alias{remove_gaps}

\title{Extract and translate codons from a set of DNA sequences}

\description{
  The \code{translate_codons()} function extracts and translates codons
  from a set of DNA sequences.
}

\usage{
translate_codons(dna, offset=0, with.init.codon=FALSE)

## Used internally by translate_codons():
extract_codons(dna, offset=0)
remove_gaps(dna, gap_letter=".")
}

\arguments{
  \item{dna}{
    A \link[Biostrings]{DNAStringSet} (or \link[Biostrings]{DNAString})
    object containing the codons to translate.

    Note that if the sequences in \code{dna} contain gaps (represented
    by the \code{.} letter) then \code{translate_codons()} and
    \code{extract_codons()} will remove them (by calling
    \code{remove_gaps()} internally) before doing anything else.
  }
  \item{offset}{
    The number of nucleotides that precede the first codon to translate.
    This must be supplied as a numeric vector with one value per
    sequence in \code{dna}, or as a single value. If the latter,
    then the same offset is used for all sequences.
  }
  \item{with.init.codon}{
    Is the first codon to translate in each DNA sequence the initiation
    codon? By default, \code{with.init.codon} is set to \code{FALSE}, in
    which case \code{translate_codons()} assumes that the first codon
    to translate in each DNA sequence is \emph{not} the initiation codon.

    See documentation of the \code{no.init.codon} argument in
    \code{?\link[Biostrings]{translate}} in the \pkg{Biostrings} package
    for more information.
  }
  \item{gap_letter}{
    The letter or symbol representing gaps. Note that IMGT and
    AIRR-community/OGRDB typically use dots (\code{"."}) to represent
    gaps in germline V gene allele sequences.
  }
}

\value{
  \code{translate_codons()} returns:
  \itemize{
    \item An \link[Biostrings]{AAStringSet} object with one amino acid
          sequence per input sequence if a \link[Biostrings]{DNAStringSet}
          object was supplied.
    \item An \link[Biostrings]{AAString} object if a
          \link[Biostrings]{DNAString} object was supplied.
  }

  \code{extract_codons()} returns:
  \itemize{
    \item A \link[Biostrings]{DNAStringSet} object with one sequence per
          input sequence if a \link[Biostrings]{DNAStringSet} object was
          supplied.
    \item A \link[Biostrings]{DNAString} object if a
          \link[Biostrings]{DNAString} object was supplied.
  }
  The sequences returned by \code{extract_codons()} are obtained by
  (1) removing all gaps (i.e. all \code{.} letters) from the input sequences,
  and (2) trimming the gap-free sequences as follow:
  \itemize{
    \item On their 5' end, sequences are trimmed by the amount of
          nucleotides specified in \code{offset}.
    \item On their 3' end, sequences are trimmed by the smallest
          amount of nucleotides that makes the length of the trimmed
          sequence a multiple of 3. Note that this will always be 0, 1,
          or 2 nucleotides.
  }

  \code{remove_gaps()} returns the gap-free versions of the input sequences
  in an object of the same type as the input object.
}

\seealso{
  \itemize{
    \item \link[Biostrings]{DNAStringSet} and \link[Biostrings]{AAStringSet}
          objects in the \pkg{Biostrings} package.

    \item The \code{\link[Biostrings]{translate}} function in the
          \pkg{Biostrings} package upon which \code{translate_codons()}
          is based.

    \item \code{\link{list_germline_dbs}} to list all the \emph{cached
          germline dbs}, that is, all the germline databases currently
          installed in \pkg{igblastr}'s persistent cache.
  }
}

\examples{
## ---------------------------------------------------------------------
## translate_codons()
## ---------------------------------------------------------------------

## Gaps will be ignored:
dna1 <- DNAStringSet(c(".TTGTCCTTTAT..A", "GAAT.CATTTATC", "CTGTCGTTTATT"))
translate_codons(dna1)

## Handling of initiation codons (see '?translate' in the Biostrings
## package for how initiation codons are handled):
translate_codons(dna1, with.init.codon=TRUE)

## Load germline V gene allele sequences for human:
list_germline_dbs()
db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
V_alleles <- load_germline_db(db_name, region_types="V")
V_alleles  # DNAStringSet object

## Translate them:
V_aa <- translate_codons(V_alleles)
V_aa  # AAStringSet object

## Some human germline V gene allele sequences have a stop codon:
has_stop_codon <- grepl("*", as.character(V_aa), fixed=TRUE)
V_aa[has_stop_codon]

## ---------------------------------------------------------------------
## extract_codons()
## ---------------------------------------------------------------------
extract_codons(dna1)
extract_codons(dna1, offset=1)

dna2 <- DNAStringSet(c("CCCAAAGGGTTT",
                       "CCAAAGGGTTT",
                       "CAAAGGGTTT",
                       "AAAGGGTTT"))
extract_codons(dna2)
extract_codons(dna2, offset=1)
extract_codons(dna2, offset=2)
extract_codons(dna2, offset=3)
extract_codons(dna2, offset=4)

extract_codons(dna2, offset=3:0)
extract_codons(dna2, offset=4:1)
extract_codons(dna2, offset=5:2)
extract_codons(dna2, offset=6:3)

## ---------------------------------------------------------------------
## remove_gaps()
## ---------------------------------------------------------------------
remove_gaps(dna1)
}

\keyword{utilities}

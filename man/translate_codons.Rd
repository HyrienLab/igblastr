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
remove_gaps(sequences, gap_letter=".")
}

\arguments{
  \item{dna}{
    A \link[Biostrings]{DNAStringSet} (or \link[Biostrings]{DNAString})
    object containing the codons to translate.
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
  \item{sequences}{
    A \link[Biostrings]{DNAStringSet} or \link[Biostrings]{AAStringSet}
    object, or other \link[Biostrings]{XStringSet} derivative.
    Can also be a \link[Biostrings]{DNAString} or \link[Biostrings]{AAString}
    object, or other \link[Biostrings]{XString} derivative.
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
  The sequences returned by \code{extract_codons()} are obtained by trimming
  the input sequences as follow:
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

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.
  }
}

\examples{
## ---------------------------------------------------------------------
## translate_codons(): basic examples
## ---------------------------------------------------------------------

## Handling of gaps (see '?translate' in the Biostrings package for how
## gaps are handled):
dna1 <- DNAString("..........TTACTCAGCCAAGCTCT...GTGTCTACGTCTCTAGGA")
translate_codons(dna1)

## Handling of initiation codons (see '?translate' in the Biostrings
## package for how initiation codons are handled):
dna2 <- DNAString("TTGATATGGCCC......TTATAA")
translate_codons(dna2)
translate_codons(dna2, with.init.codon=TRUE)

## ---------------------------------------------------------------------
## translate_codons(): advanced example
## ---------------------------------------------------------------------
if (IMGT_is_up()) {
  dir.create(destdir <- tempfile("fasta_dir_"))
  download_IMGT_germline_sequences("202631-1", organism="Homo sapiens",
                                   destdir=destdir)
  gapped_V_alleles <- readDNAStringSet(file.path(destdir, "IGKV.fasta"))

  ## Translate the gapped V alleles:
  gapped_V_alleles_aa <- translate_codons(gapped_V_alleles)
  gapped_V_alleles_aa

  ## Cysteins are expected at positions 23 and 104 (see
  ## https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html):
  table(substr(gapped_V_alleles_aa, 23, 23))
  table(substr(gapped_V_alleles_aa, 104, 104))

  unlink(destdir, recursive=TRUE)
}

## ---------------------------------------------------------------------
## extract_codons()
## ---------------------------------------------------------------------
dna3 <- DNAStringSet(c("CCCAAAGGGTTT",
                       "CCAAAGGGTTT",
                       "CAAAGGGTTT",
                       "AAAGGGTTT"))
extract_codons(dna3)
extract_codons(dna3, offset=1)
extract_codons(dna3, offset=2)
extract_codons(dna3, offset=3)
extract_codons(dna3, offset=4)

extract_codons(dna3, offset=3:0)
extract_codons(dna3, offset=4:1)
extract_codons(dna3, offset=5:2)
extract_codons(dna3, offset=6:3)

## ---------------------------------------------------------------------
## remove_gaps()
## ---------------------------------------------------------------------
remove_gaps(dna1)
}

\keyword{utilities}

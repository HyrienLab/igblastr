\name{compute_auxdata}

\alias{compute_auxdata}

\title{Compute IgBLAST auxiliary data}

\description{
  A utility function to annotate a set of germline J gene allele sequences
  in a way similar to how they are annotated in IgBLAST \emph{auxiliary data}.

  Note that the annotation produced by the function can be used by IgBLAST
  as a substitute to the \emph{auxiliary data} shipped with the IgBLAST
  software (and typically included in a standard IgBLAST installation).

  See \code{?\link{load_auxdata}} for more information about IgBLAST
  \emph{auxiliary data}.
}

\usage{
compute_auxdata(J_alleles, codon_starts=NULL, no.warning=FALSE)
}

\arguments{
  \item{J_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline J
    gene allele sequences.
  }
  \item{codon_starts}{
    When supplied, must be a named integer vector \emph{parallel} to
    \code{J_alleles}, that is, it must have the same length and names
    as \code{J_alleles}. For each allele in \code{J_alleles},
    \code{codon_starts} should indicate the start position (1-based) of
    the first codon in the allele. This should be a number >= 1 and <= 3.
    \code{NA}s are allowed when the coding frame is not known.

    Note that the "coding frame start" is 0-based so can simply be inferred
    from the "codon start" by subtracting 1.

    \code{compute_auxdata()} will use the supplied "codon starts" to
    either validate the "coding frame starts" that it could determine,
    or to set the "coding frame starts" that it could not determine.
  }
  \item{no.warning}{
    Set to \code{TRUE} to suppress the warning that \code{compute_auxdata()}
    otherwise issues when the \code{coding_frame_start} column of the
    returned data.frame contains \code{NA}s.
  }
}

\details{
  The FWR4 region is expected to start with the following amino acid
  motifs (X represents any amino acid):
  \itemize{
    \item "WGXG" on the heavy chain;
    \item "FGXG" on the light chain.
  }

  \code{compute_auxdata()} searches for the "WGXG" and "FGXG" motifs in
  the supplied allele sequences to determine the start of their FWR4
  region. From there it can easily infer the \code{cdr3_end},
  \code{coding_frame_start}, and \code{extra_bps} columns.

  Note that the function will emit a warning if the start of the FWR4
  region (and therefore the CDR3 end) could not be found for some alleles,
  or if a stop codon was found in some alleles.
}

\value{
  Returns the computed \emph{auxiliary data} in a data.frame with
  1 row per supplied J allele sequence and the same columns as the
  data.frame returned by \code{read_auxdata()}.
  See \code{?\link{read_auxdata}} for more information.
}

\seealso{
  \itemize{
    \item \code{\link{load_auxdata}} to access IgBLAST \emph{auxiliary data}.

    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.

    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about IgBLAST \emph{auxiliary data}.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------

## Let's load a set of human J allele sequences:
db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
J_alleles <- load_germline_db(db_name, region_types="J")
J_alleles  # DNAStringSet object

computed_auxdata <- compute_auxdata(J_alleles)
head(computed_auxdata)

## ---------------------------------------------------------------------
## SANITY CHECK
## ---------------------------------------------------------------------

## Note that 'computed_auxdata' is in agreement with the auxiliary
## data included in IgBLAST for human, except for alleles IGHJ6*02
## and IGHJ6*03 (the 'extra_bps' column contains incorrect values
## for these two alleles, 0's instead of 1's):
human_auxdata0 <- load_auxdata("human", which="original")
bad_alleles <- c("IGHJ6*02", "IGHJ6*03")
subset(human_auxdata0, allele_name \%in\% bad_alleles)

## We manually correct this:
fixme <- human_auxdata0[ , "allele_name"] \%in\% bad_alleles
human_auxdata0[fixme, "extra_bps"] <- 1L  # replace 0L with 1L

## Now the data in 'computed_auxdata' matches exactly the corresponding
## data in 'human_auxdata0':
m <- match(computed_auxdata[ , "allele_name"],
           human_auxdata0[ , "allele_name"])
human_auxdata <- human_auxdata0[m, ]
rownames(human_auxdata) <- NULL
stopifnot(identical(computed_auxdata, human_auxdata))
}

\keyword{utilities}

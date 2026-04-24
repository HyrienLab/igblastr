\name{auxdata-IO}

\alias{auxdata-IO}
\alias{auxdata_IO}

\alias{read_auxdata}
\alias{write_auxdata}

\title{Read/write IgBLAST auxiliary data files}

\description{
  IgBLAST auxiliary data files (a.k.a. .aux files) are tab-delimited text
  files used by IgBLAST to annotate germline J gene allele sequences.

  The \pkg{igblastr} package provides low-level functions to read/write
  this type of file.
}

\usage{
read_auxdata(filepath)
write_auxdata(auxdata, file="")
}

\arguments{
  \item{filepath}{
    The path to the IgBLAST auxiliary data file to read.
  }
  \item{auxdata}{
    A data.frame with 1 row per germline J gene allele sequence and the
    same columns as the data.frame returned by \code{read_auxdata()}.
    See Value section below for the details.
  }
  \item{file}{
    The path to the IgBLAST auxiliary data file to write.
  }
}

\value{
  \code{read_auxdata()} returns a data.frame with 1 row per germline J gene
  allele sequence and the following columns:
  \enumerate{
    \item \code{allele_name}: allele name;
    \item \code{coding_frame_start}: first coding frame start
          position (0-based);
    \item \code{chain_type}: chain type as a 2-letter code e.g. JK
          for "J allele from the Kappa locus" (BCR locus) or JG
          for "J allele from the Gamma locus" (TCR locus);
    \item \code{cdr3_end}: CDR3 end position (0-based);
    \item \code{extra_bps}: extra base pairs beyond J coding end.
  }

  \code{write_auxdata()} doesn't return anything (i.e. invisible \code{NULL}).
}

\seealso{
  \itemize{
    \item \code{\link{compute_auxdata}} to annotate a set of
          germline J gene allele sequences.

    \item \code{\link{extract_auxdata_from_ogrdb_json}} to extract the
          coding frame and CDR3 end information for the J alleles annotated
          in an AIRR-C JSON file.

    \item \link{auxdata_utils} to access IgBLAST \emph{auxiliary data}.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
## COMING SOON...
}

\keyword{utilities}

\name{ndm_data-IO}

\alias{ndm_data-IO}
\alias{ndm_data_IO}

\alias{read_ndm_data}
\alias{write_ndm_data}
\alias{validate_ndm_rows}
\alias{same_ndm_data}

\title{Read/write "ndm" files}

\description{
  The \emph{ndm} (Nucleotide Domain Mapping) file format is a specialized,
  tab-delimited text file used by IgBLAST to define the FWR/CDR boundaries
  on the germline V gene allele sequences.

  The \pkg{igblastr} package provides low-level functions to read/write data
  in \emph{ndm} format.
}

\usage{
read_ndm_data(filepath)
write_ndm_data(ndm_data, file="", check.data=FALSE)

validate_ndm_rows(ndm_data, allow.repeated.rows=FALSE)
}

\arguments{
  \item{filepath}{
    The path to the "ndm" file to read.
  }
  \item{ndm_data}{
    A data.frame with 1 row per germline V gene allele sequence and the
    same columns as the data.frame returned by \code{read_ndm_data()}.
    See Value section below for the details.
  }
  \item{file}{
    The path to the "ndm" file to write.
  }
  \item{check.data}{
    COMING SOON...
  }
  \item{allow.repeated.rows}{
    COMING SOON...
  }
}

\value{
  \code{read_ndm_data()} returns a data.frame with 1 row per germline V gene
  allele sequence and the following columns:
  \itemize{
    \item \code{allele_name}: allele name;
    \item \code{fwr1_start}, \code{fwr1_end}: FWR1 start/end
          positions (1-based);
    \item \code{cdr1_start}, \code{cdr1_end}: CDR1 start/end
          positions (1-based);
    \item \code{fwr2_start}, \code{fwr2_end}: FWR2 start/end
          positions (1-based);
    \item \code{cdr2_start}, \code{cdr2_end}: CDR2 start/end
          positions (1-based);
    \item \code{fwr3_start}, \code{fwr3_end}: FWR3 start/end
          positions (1-based);
    \item \code{chain_type}: chain type as a 2-letter code e.g. VK
          for "V allele from the Kappa locus" (BCR locus) or VG
          for "V allele from the Gamma locus" (TCR locus);
    \item \code{coding_frame_start}: first coding frame start
          position (0-based).
  }

  \code{write_ndm_data()} doesn't return anything (i.e. invisible \code{NULL}).

  \code{validate_ndm_rows()} returns a logical vector with 1 value per row
  in the \code{ndm_data} data.frame indicating whether the data in that
  row is valid or not.
}

\seealso{
  \itemize{
    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.

    \item \code{\link{makeogrannote}} to extract the germline V gene allele
          FWR/CDR annotations (a.k.a. V gene delineations) from an
          AIRR-C JSON file.

    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

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

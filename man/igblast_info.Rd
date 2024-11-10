\name{igblast_info}

\alias{has_igblast}
\alias{igblast_info}
\alias{igblast_version}

\alias{print.igblast_info}

\title{Check IgBLAST used by igblastr}

\description{
  Collect basic information about the IgBLAST installation used by
  the \pkg{igblastr} package.
}

\usage{
has_igblast()
igblast_info()
igblast_version()
}

\value{
  \code{has_igblast()} returns \code{TRUE} or \code{FALSE}.

  \code{igblast_info()} returns a named list containing various
  information about the IgBLAST used by \pkg{igblastr}.

  \code{igblast_version()} returns a single string containing the
  version of the IgBLAST used by \pkg{igblastr}.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          program from R. This is the main function in the \pkg{igblastr}
          package.

    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \code{\link{get_igblast_root}} to get (or set) the IgBLAST
          installation used (or to be used) by the \pkg{igblastr} package.

    \item \link{IGBLAST_ROOT} to set the \emph{external} IgBLAST installation
          to be used by the \pkg{igblastr} package in a persistent manner.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()
}

\keyword{utilities}

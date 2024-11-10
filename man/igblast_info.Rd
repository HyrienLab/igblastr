\name{igblast_info}

\alias{igblast_info}

\alias{has_igblast}
\alias{igblast_version}
\alias{list_igblast_organisms}

\alias{print.igblast_info}

\title{Check IgBLAST used by igblastr}

\description{
  Collect basic information about the IgBLAST installation used by
  the \pkg{igblastr} package.
}

\usage{
igblast_info()

has_igblast()
igblast_version()
list_igblast_organisms()
}

\value{
  \code{igblast_info()} returns a named list containing basic
  information about the IgBLAST installation used by \pkg{igblastr}.

  \code{has_igblast()} returns \code{TRUE} or \code{FALSE}.

  \code{igblast_version()} returns a single string containing the
  version of the IgBLAST used by \pkg{igblastr}.

  \code{list_igblast_organisms()} returns a character vector that
  lists the organisms for which IgBLAST provides internal data.
  Note that this is obtained by simply listing the content of the
  \code{internal_data} directory in the IgBLAST installation
  used by \pkg{igblastr}.
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

list_igblast_organisms()
}

\keyword{utilities}

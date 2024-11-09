\name{igblast_info}

\alias{igblast_info}
\alias{igblast_version}
\alias{get_igblast_root}
\alias{set_igblast_root}

\alias{print.igblast_info}

\title{Check IgBLAST used by igblastr}

\description{
  Get and report information about the IgBLAST installation used
  by the \pkg{igblastr} package.
}

\usage{
igblast_info()
igblast_version()

get_igblast_root()
set_igblast_root(version_or_path)
}

\arguments{
  \item{version_or_path}{
    A single string that contains either a version number (e.g.
    \code{"1.22.0"}) or the path to an IgBLAST installation.
  }
}

\details{
  \code{set_igblast_root} can be used to set or change the path to
  the IgBLAST installation to use. This can be an \emph{internal} or
  \emph{external} installation.

  In the former case, \code{version_or_path} should be the version
  of an existing \emph{internal} installation. The setting will be
  persistent.

  In the latter case, it should be the full path (absolute or relative)
  to the \emph{root directory} of a valid \emph{external} installation.
  Note that the setting won't be persistent i.e. it won't be
  remembered across R sessions. See \code{?IGBLAST_ROOT} for how to set
  the \emph{external} IgBLAST installation to use in \pkg{igblastr} in
  a persistent manner.
}

\value{
  \code{igblast_info()} returns a named list containing various
  information about the IgBLAST used by \pkg{igblastr}.

  \code{igblast_version()} returns a single string containing the
  version of the IgBLAST used by \pkg{igblastr}.

  \code{get_igblast_root()} returns a single string containing the path
  to the \emph{root directory} of the IgBLAST installation used by
  \pkg{igblastr}.

  \code{set_igblast_root()} returns a single string containing the path
  to the \emph{root directory} of the newly selected IgBLAST installation.
  The string is returned invisibly.
}

\seealso{
  \itemize{
    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \link{IGBLAST_ROOT} to select the \emph{external} IgBLAST
          installation to use in \pkg{igblastr} in a persistent manner.

    \item The \code{\link{igblastn}} function to invoke the igblastn
          command from R.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
igblast_info()

get_igblast_root()
}

\keyword{utilities}

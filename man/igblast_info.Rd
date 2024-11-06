\name{igblast_info}

\alias{igblast_info}
\alias{igblast_version}
\alias{get_igblast_root}
\alias{set_igblast_root}

\alias{print.igblast_info}

\title{Check IgBLAST used by igblastr}

\description{
  Get and report information about the IgBLAST installation used
  by \pkg{igblastr}.
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
  \code{set_igblast_root} can be used to set the path to the IgBLAST
  installation to use. This can be either an \emph{internal} or
  an \emph{external} installation.

  In the former case, \code{version_or_path} should be the version
  of an existing \emph{internal} installation. The setting will be
  persistent.

  In the latter case, it should be the full path (absolute or relative)
  to the \emph{root directory} of a valid \emph{external} installation.
  Note that the setting won't be persistent, that is, it won't be
  remembered when you quit R. See \code{?IGBLAST_ROOT} for how to set
  the \emph{external} IgBLAST installation to use in \pkg{igblastr} in
  a persistent manner.
}

\value{
  \code{igblast_info()} returns the information about the IgBLAST used
  by \pkg{igblastr} in a named list.

  \code{igblast_version()} returns the version of the IgBLAST used
  by \pkg{igblastr} in a single string.

  \code{get_igblast_root()} returns the path to the \emph{root directory}
  of the IgBLAST installation used by \pkg{igblastr}.

  \code{set_igblast_root()} returns the path to the \emph{root directory}
  of the newly selected IgBLAST installation, as an invisible string.
}

\seealso{
  \itemize{
    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \link{IGBLAST_ROOT} to select the \emph{external} IgBLAST
          installation to use in \pkg{igblastr} in a persistent manner.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
igblast_info()

get_igblast_root()
}

\keyword{utilities}

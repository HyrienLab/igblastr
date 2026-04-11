\name{igblast_info}

\alias{igblast_info}

\alias{igblast_build}
\alias{igblastn_version}
\alias{makeblastdb_version}
\alias{list_igblast_organisms}
\alias{has_igblast}

\alias{print.igblast_info}

\title{Check IgBLAST used by igblastr}

\description{
  Collect basic information about the IgBLAST installation used by
  the \pkg{igblastr} package, or about any IgBLAST installation on
  the user machine.
}

\usage{
igblast_info(igblast_root=get_igblast_root())

igblast_build(igblast_root=get_igblast_root())
igblastn_version(igblast_root=get_igblast_root(), raw.version=FALSE)
makeblastdb_version(igblast_root=get_igblast_root(), raw.version=FALSE)
list_igblast_organisms(igblast_root=get_igblast_root())

has_igblast()
}

\arguments{
  \item{igblast_root}{
    A single string that is the path to an IgBLAST installation.
    By default \code{igblast_root} is set to \code{get_igblast_root()},
    which is the path to the IgBLAST installation used by the \pkg{igblastr}
    package. See \code{?\link{get_igblast_root}} for more information.

    Note that the supplied string must contain the path to the
    \emph{root directory} of an IgBLAST installation, that is, to a directory
    with a \code{bin} subdirectory in it that has the \code{igblastn},
    \code{igblastp}, and \code{makeblastdb} \emph{standalone executables}
    (on Windows these executables are files named \code{igblastn.exe},
    \code{igblastp.exe}, and \code{makeblastdb.exe}, respectively).
  }
  \item{raw.version}{
    By default (i.e. when \code{raw.version} is omitted or set to
    \code{FALSE}), \code{igblastn_version()} and \code{makeblastdb_version()}
    return the version string of the \code{igblastn} and \code{makeblastdb}
    \emph{standalone executables} included in IgBLAST. This string is
    extracted from the output produced by system commands:
    \preformatted{    igblastn -version}
    and
    \preformatted{    makeblastdb -version}

    When \code{raw.version} is set to \code{TRUE}, \code{igblastn_version()}
    and \code{makeblastdb_version()} return the \emph{full ouput} produced
    by the above commands.
  }
}

\value{
  \code{igblast_info()} returns a named list containing basic
  information about the IgBLAST installation.

  \code{igblast_build()} returns a single string containing IgBLAST build
  information.

  By default, \code{igblastn_version()} returns a single string containing
  the version of the \code{igblastn} \emph{standalone executable} included
  in IgBLAST.

  By default, \code{makeblastdb_version()} returns a single string containing
  the version of the \code{makeblastdb} \emph{standalone executable} included
  in IgBLAST.

  \code{list_igblast_organisms()} returns a character vector that
  lists the organisms for which IgBLAST provides \emph{internal data}.
  Note that this is obtained by simply listing the content of the
  \code{internal_data} directory located in the IgBLAST installation.

  \code{has_igblast()} returns \code{TRUE} or \code{FALSE}.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \code{\link{get_igblast_root}} to get (or set) the IgBLAST
          installation used (or to be used) by the \pkg{igblastr} package.

    \item \link{IGBLAST_ROOT} to set the \emph{external} IgBLAST installation
          to be used by the \pkg{igblastr} package in a persistent manner.

    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

    \item \link{auxdata_utils} to access IgBLAST \emph{auxiliary data}.

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

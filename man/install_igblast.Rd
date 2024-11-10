\name{install_igblast}

\alias{install_igblast}

\title{Install IgBLAST}

\description{
  Download and install a pre-compiled IgBLAST from NCBI FTP site for
  use with \pkg{igblastr}.
}

\usage{
install_igblast(release="LATEST", force=FALSE, ...)
}

\arguments{
  \item{release}{
    A single string specifying the IgBLAST release version to install.
    For example \code{"LATEST"} (recommended), or one
    of the IgBLAST release versions listed at
    \url{https://ftp.ncbi.nih.gov/blast/executables/igblast/release/}
    (e.g. \code{"1.21.0"}).
    Note that old versions have not been tested and are not
    guaranteed to be compatible with the \pkg{igblastr} package.
  }
  \item{force}{
    Set to \code{TRUE} to reinstall if the specified IgBLAST release version
    is already installed.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?\link[utils]{download.file}} in the \pkg{utils} package for
    more information.
  }
}

\value{
  The path to the \emph{root directory} of the IgBLAST installation,
  as an invisible string.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          program from R. This is the main function in the \pkg{igblastr}
          package.

    \item \link{IGBLAST_ROOT} to use an \emph{external} IgBLAST installation.

    \item \code{\link{igblast_info}} to collect basic information about
          the IgBLAST installation used by the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()
}

\keyword{utilities}

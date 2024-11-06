\name{install_igblast}

\alias{install_igblast}

\title{Install IgBlast}

\description{
  Download and install a pre-compiled IgBlast from NCBI FTP site for
  use with \pkg{igblastr}.
}

\usage{
install_igblast(release="LATEST", force=FALSE, ...)
}

\arguments{
  \item{release}{
    A single string specifying the IgBlast release version to install.
    For example \code{"LATEST"} (recommended), or one
    of the IgBlast release versions listed at
    \url{ftp.ncbi.nih.gov/blast/executables/igblast/release/}
    (e.g. \code{"1.21.0"}).
    Note that old versions have not been tested and are not
    guaranteed to be compatible with the \pkg{igblastr} package.
  }
  \item{force}{
    Set to \code{TRUE} to reinstall if the specified IgBlast release version
    is already installed.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?download.file} in the \pkg{utils} package for more information.
  }
}

\value{
  The path to the root of the IgBlast installation, as an invisible string.
}

\seealso{
  \itemize{
    \item \link{IGBLAST_ROOT} to use an \emph{external} IgBlast installation.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
install_igblast(force=TRUE)

igblast_info()
}

\keyword{utilities}

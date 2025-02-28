\name{auxiliary-data-utils}

\alias{auxiliary-data-utils}
\alias{auxiliary_data_utils}

\alias{get_igblast_auxiliary_data}

\title{Manipulation of IgBLAST auxiliary data}

\description{
  A standard IgBLAST installation -- like the one used by the \pkg{igblastr}
  package -- typically includes \emph{auxiliary data} for various
  organisms, in the form of one tabulated file per organism. Each file
  indicates the germline J gene coding frame start position, the J gene type,
  and the CDR3 end position for each sequence in the germline J sequence
  database. See \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
  for the details.

  You can use \code{get_igblast_auxiliary_data()} to obtain the path to the
  file containing the auxiliary data for a given organism.
}

\usage{
get_igblast_auxiliary_data(organism)
}

\arguments{
  \item{organism}{
    A single string containing the name of an organism as
    returned by \code{\link{list_igblast_organisms}()}.
  }
}

\details{
  WARNING: According to
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
  the auxiliary data included in IgBLAST is specific to a particular
  NCBI or IMGT germline db. Unfortunately this means that this data is
  NOT guaranteed to be compatible with the germline db that you will
  use with \code{\link{igblastn}()}. See documentation of the
  \code{auxiliary_data} argument in \code{?\link{igblastn}} for
  more information about this.
}

\value{
  Returns a single string containing the path to the auxiliary data
  included in the IgBLAST installation used by \pkg{igblastr}, for
  the specified organism. Not necessarily suitable to use with
  \code{\link{igblastn}()} (see above WARNING).
}

\seealso{
  \itemize{
    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about the IgBLAST auxiliary data.

    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          executable from R. This is the main function in the \pkg{igblastr}
          package.

    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \code{\link{get_igblast_root}} to get (or set) the IgBLAST
          installation used (or to be used) by the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

list_igblast_organisms()

get_igblast_auxiliary_data("human")

## Make sure to read the WARNING above before using this auxiliary
## data with igblastn().
}

\keyword{utilities}

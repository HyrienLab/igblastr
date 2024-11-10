\name{igblastn}

\alias{igblastn}
\alias{igblastn_help}

\title{Run the igblastn program from R}

\description{
  \code{igblastn()} is the main function in the \pkg{igblastr} package.
}

\usage{
igblastn(query, outfmt="AIRR", organism="auto", ...,
         show.in.browser=FALSE, show.command.only=FALSE)

igblastn_help(long.help=FALSE, show.in.browser=FALSE)
}

\arguments{
  \item{query}{
    TODO
  }
  \item{outfmt}{
    TODO
  }
  \item{organism}{
    TODO
  }
  \item{...}{
    Extra arguments to be passed to the igblastn program.
  }
  \item{show.in.browser}{
    TODO
  }
  \item{show.command.only}{
    TODO
  }
  \item{long.help}{
    TODO
  }
}

\value{
  TODO
}

\seealso{
  \itemize{
    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \code{\link{igblast_info}} to collect basic information about
          the IgBLAST installation used by the \pkg{igblastr} package.

    \item \code{\link{install_VQUEST_germline_db}} to install a V-QUEST
          germline db.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## File 'catnap_bnabs.fasta' contains a small BCR sequencing data set
## made of 1000 heavy- and light-chains associated with bnAbs downloaded
## from the CATNAP database:
catnap_bnabs <- system.file(package="igblastr", "extdata",
                            "catnap_bnabs.fasta")

## Install germline db from V-QUEST:
db_name <- install_VQUEST_germline_db("202430-2", "Homo_sapiens",
                                      db_type="IG", force=TRUE)

## Select germline db to use with igblastn():
use_germline_db(db_name)

## Run igblastn:
igblastn(catnap_bnabs, show.in.browser=TRUE)
}

\keyword{manip}

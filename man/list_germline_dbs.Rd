\name{list_germline_dbs}

\alias{list_germline_dbs}
\alias{use_germline_db}

\title{List local germline dbs and select one to use with igblastn()}

\description{
  Use \code{list_germline_dbs()} to list all the \emph{local germline dbs},
  that is, all the germline databases currently installed in \pkg{igblastr}'s
  persistent cache.

  Use \code{use_germline_db()} to select the local germline db to use with
  \code{\link{igblastn}()}. This choice will be remembered across R sessions.
}

\usage{
list_germline_dbs(names.only=FALSE)

use_germline_db(db_name=NULL)
}

\arguments{
  \item{names.only}{
    By default \code{list_germline_dbs()} returns the list of local
    germline dbs in a data.frame with one db per row. Set \code{names.only}
    to \code{TRUE} to return only the db names in a character vector.
  }
  \item{db_name}{
    \code{NULL} or a single string specifying the name of the local germline
    db to use. Use \code{list_germline_dbs()} to list all the local germline
    dbs.

    If set to \code{NULL} (the default), then \code{use_germline_db()}
    returns the name of the local germline db that is currently in use,
    if any. Otherwise it raises an error.
  }
}

\details{
  The \pkg{igblastr} package provides a small set of utilities to
  manage the local germline and C-region databases to use with
  \code{\link{igblastn}()}.

  A \emph{local germline db} contains the nucleotide sequences of the
  V-, D-, and J-regions for a given organism.

  A \emph{local C-region db} contains the nucleotide sequences of the
  C-regions for a given organism.

  This man page documents the two core utilities to manage local germline dbs:
  \code{list_germline_dbs()} and \code{use_germline_db()}.

  The core utilities to manage local C-region dbs are documented in the man
  page for \code{\link{list_c_region_dbs}}.
}

\value{
  \code{list_germline_dbs()} returns the list of all local germline dbs in a
  data.frame with one db per row (if \code{names.only} is \code{FALSE}, the
  default), or in a character vector (if \code{names.only} is \code{TRUE}).

  When called with no argument, \code{use_germline_db()} returns a single
  string containing the name of the local germline db currently used by
  \code{\link{igblastn}()} if any, or it raises an error if no germline db
  has been selected yet.

  When called with the \code{db_name} argument, \code{use_germline_db(db_name)}
  returns \code{db_name} invisibly.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          program from R. This is the main function in the \pkg{igblastr}
          package.

    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item \code{\link{use_c_region_db}} to select the local C-region db
          to use with \code{igblastn()}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## Install Mouse germline db from IMGT (Perl required!):
install_IMGT_germline_db("202430-2", "Mus_musculus", force=TRUE)

list_germline_dbs()

list_germline_dbs(TRUE)  # return only the db names

use_germline_db("IMGT-202430-2.Mus_musculus.IG")

use_germline_db()
}

\keyword{utilities}

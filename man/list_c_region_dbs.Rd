\name{list_c_region_dbs}

\alias{list_c_region_dbs}
\alias{use_c_region_db}

\title{List local C-region dbs and select one to use with igblastn()}

\description{
  Use \code{list_c_region_dbs()} to list all the \emph{local C-region dbs},
  that is, all the C-region databases currently installed in \pkg{igblastr}'s
  persistent cache.

  Use \code{use_c_region_db()} to select the local C-region db to use with
  \code{\link{igblastn}()}. This choice will be remembered across R sessions.
}

\usage{
list_c_region_dbs(names.only=FALSE)

use_c_region_db(db_name=NULL)
}

\arguments{
  \item{names.only}{
    By default \code{list_c_region_dbs()} returns the list of local
    C-region dbs in a data.frame with one db per row. Set \code{names.only}
    to \code{TRUE} to return only the db names in a character vector.
  }
  \item{db_name}{
    \code{NULL} or a single string specifying the name of the local C-region
    db to use. Use \code{list_c_region_dbs()} to list all the local C-region
    dbs.

    If set to \code{NULL} (the default), then \code{use_c_region_db()}
    returns the name of the local C-region db that is currently in use,
    if any. Otherwise it returns the empty string (\code{""}).

    Note that the current selection can be cancelled with
    \code{use_c_region_db("")}.
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

  This man page documents the two core utilities to manage local C-region dbs:
  \code{list_c_region_dbs()} and \code{use_c_region_db()}.

  The core utilities to manage local germline dbs are documented in the man
  page for \code{\link{list_germline_dbs}}.
}

\value{
  \code{list_c_region_dbs()} returns the list of all local C-regions dbs in a
  data.frame with one db per row (if \code{names.only} is \code{FALSE}, the
  default), or in a character vector (if \code{names.only} is \code{TRUE}).

  When called with no argument, \code{use_c_region_db()} returns a single
  string containing the name of the local C-region db currently used by
  \code{\link{igblastn}()} if any, or the empty string (\code{""}) if
  \code{\link{igblastn}()} is not using any C-region db.

  When called with the \code{db_name} argument, \code{use_c_region_db(db_name)}
  returns \code{db_name} invisibly.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          program from R. This is the main function in the \pkg{igblastr}
          package.

    \item \code{\link{use_germline_db}} to select the local germline db
          to use with \code{igblastn()}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

list_c_region_dbs()

list_c_region_dbs(TRUE)  # return only the db names

use_c_region_db("IMGT.human.IG")

use_c_region_db()

use_c_region_db("")  # cancel the current selection

use_c_region_db()
}

\keyword{utilities}

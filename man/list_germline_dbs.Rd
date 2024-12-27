\name{list_germline_dbs}

\alias{list_germline_dbs}
\alias{use_germline_db}

\title{Select germline db to use with igblastn()}

\description{
  Use \code{list_germline_dbs()} to list all the germline databases that
  are installed in the cache (a.k.a. local germline dbs).

  Use \code{use_germline_db()} to select the local germline db to use with
  \code{\link{igblastn}()}. This choice will be remembered across R sessions.
}

\usage{
list_germline_dbs()

use_germline_db(db_name=NULL)
}

\arguments{
  \item{db_name}{
    \code{NULL} or a single string specifying the name of the local germline
    db to use. Use \code{list_germline_dbs()} to list all the local germline
    dbs.

    If set to \code{NULL} (the default), then \code{use_germline_db()}
    returns the name of the local germline db that is currently in use,
    if any. Otherwise it raises an error.
  }
}

\value{
  \code{list_germline_dbs()} returns the names of all local germline dbs in
  a character vector.

  \code{use_germline_db()} returns a single string containing the name of
  the local germline db in use with \code{\link{igblastn}()}.
  If called with no argument, this is the name of the db that is currently
  in use. If \code{db_name} is supplied, this is \code{db_name}. However note
  that in this case the string is returned invisibly.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          program from R. This is the main function in the \pkg{igblastr}
          package.

    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## Install Mouse germline db from IMGT (Perl required!):
install_IMGT_germline_db("202430-2", "Mus_musculus", force=TRUE)

list_germline_dbs()

use_germline_db("IMGT-202430-2.Mus_musculus.IG")

use_germline_db()
}

\keyword{utilities}

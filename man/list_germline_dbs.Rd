\name{list_germline_dbs}

\alias{list_germline_dbs}
\alias{use_germline_db}
\alias{load_germline_db}

\title{List local germline dbs and select one to use with igblastn()}

\description{
  Use \code{list_germline_dbs()} to list all the \emph{local germline dbs},
  that is, all the germline databases currently installed in \pkg{igblastr}'s
  persistent cache.

  Use \code{use_germline_db()} to select the local germline db to use with
  \code{\link{igblastn}()}. This choice will be remembered across R sessions.

  Use \code{load_germline_db()} to load the nucleotide sequences of the
  gene regions stored in a local germline db.
}

\usage{
list_germline_dbs(names.only=FALSE)

use_germline_db(db_name=NULL)

load_germline_db(db_name, region_types=NULL)
}

\arguments{
  \item{names.only}{
    By default \code{list_germline_dbs()} returns the list of local
    germline dbs in a data.frame with one db per row. Set \code{names.only}
    to \code{TRUE} to return only the db names in a character vector.
  }
  \item{db_name}{
    For \code{use_germline_db()}:

    \code{NULL} or a single string specifying the name of the
    local germline db to use. Use \code{list_germline_dbs()} to
    list all the local germline dbs.

    If set to \code{NULL} (the default), then
    \code{use_germline_db()} returns the name of the local
    germline db that is currently in use, if any. Otherwise it
    raises an error.

    For \code{load_germline_db()}:

    A single string specifying the name of the local germline db
    from which to load the V, D, and/or J regions. Use
    \code{list_germline_dbs()} to list all the local germline dbs.
  }
  \item{region_types}{
    The types of regions (V, D, and/or J) to load from the database.
    Specified as a single string (e.g. \code{"DJ"}) or as a character
    vector of single-letter elements (e.g. \code{c("D", "J")}).
    By default (i.e. when \code{region_types} is \code{NULL}), all the
    regions are returned.
  }
}

\details{
  The \pkg{igblastr} package provides a small set of utilities to
  manage the local germline and C-region databases to use with
  \code{\link{igblastn}()}.

  Terminology:
  \itemize{
    \item A \emph{local germline db} contains the nucleotide sequences of
          the V, D, and J gene regions for a given organism.
    \item A \emph{local C-region db} contains the nucleotide sequences of
          the C regions (i.e. constant gene regions) for a given organism.
  }

  This man page documents the basic utilities to operate on local
  germline dbs: \code{list_germline_dbs()}, \code{use_germline_db()},
  and \code{load_germline_db()}.

  The basic utilities to operate on local C-region dbs are documented
  in the man page for \code{\link{list_c_region_dbs}}.
}

\value{
  \code{list_germline_dbs()} returns the list of all local germline dbs
  in a data.frame with one db per row (if \code{names.only} is \code{FALSE},
  which is the default), or in a character vector (if \code{names.only}
  is \code{TRUE}). Columns IGHV, IGHD, IGHJ, IGKV, IGKJ, IGLV, and IGLJ
  in the data.frame indicate the number of sequences from each group.

  When called with no argument, \code{use_germline_db()} returns a single
  string containing the name of the local germline db currently used by
  \code{\link{igblastn}()} if any, or it raises an error if no germline db
  has been selected yet.

  When called with the \code{db_name} argument, \code{use_germline_db(db_name)}
  returns \code{db_name} invisibly.

  \code{load_germline_db()} returns the nucleotide sequences from the
  specified germline db in a named \link[Biostrings]{DNAStringSet} object.
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

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## Install Mouse germline db from IMGT (Perl required!):
install_IMGT_germline_db("202449-1", "Mus_musculus", force=TRUE)

list_germline_dbs()

list_germline_dbs(TRUE)  # return only the db names

db_name <- "IMGT-202449-1.Mus_musculus.IGH+IGK+IGL"
use_germline_db(db_name)

use_germline_db()

load_germline_db(db_name)
load_germline_db(db_name, region_types="D")
load_germline_db(db_name, region_types="DJ")
}

\keyword{utilities}

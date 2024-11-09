\name{install_VQUEST_germline_db}

\alias{install_VQUEST_germline_db}
\alias{list_VQUEST_releases}

\title{Install a V-QUEST germline db}

\description{
  Download V-QUEST germline sequences from the IMGT website, and create a
  small local database with the subset of sequences of a given organism.
  The small local database can then be used with \code{\link{igblastn}()}.
}

\usage{
install_VQUEST_germline_db(release, organism="Homo sapiens",
                           db_type=c("IG", "TR", "IG-TR"),
                           force=FALSE, ...)

list_VQUEST_releases(recache=FALSE)
}

\arguments{
  \item{release}{
    A single string specifying the V-QUEST release to get the germline
    sequences from.
    Use \code{list_VQUEST_releases()} to list all releases.
  }
  \item{organism}{
    A single string specifying the organism for which to get the germline
    sequences.
  }
  \item{db_type}{
    A single string specifying the type of germline sequences to store
    in the local database:
    \itemize{
      \item \code{"IG"} (the default) stands for \emph{Immunoglobulin};
      \item \code{"TR"} stands for \emph{T-cell receptor};
      \item \code{"IG-TR"} means that both types should go in the local
            database.
    }
  }
  \item{force}{
    Set to \code{TRUE} to reinstall if the requested database is already
    installed.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?download.file} in the \pkg{utils} package for more information.
  }
  \item{recache}{
    \code{list_VQUEST_releases()} uses a caching mechanism so that the
    list of V-QUEST releases gets downloaded only once from the IMGT
    website during an R session (note that this caching is done in memory
    so it does not persist across sessions).
    Set \code{recache} to \code{TRUE} to force a new download (and
    recaching) of the list of V-QUEST releases.
  }
}

\value{
  \code{install_VQUEST_germline_db()} returns the path to the local
  germline db as an invisible string.

  \code{list_VQUEST_releases()} returns the list of V-QUEST releases
  in a character vector. The releases are sorted from newest to oldest
  (latest release is first).
}

\seealso{
  \itemize{
    \item \code{\link{use_germline_db}} to select the local germline
          db to use with \code{igblastn()}.

    \item The IMGT website: \url{https://www.imgt.org/}.

    \item The V-QUEST download site:
          \url{https://www.imgt.org/download/V-QUEST/}.

    \item The \code{\link{igblastn}} function to invoke the igblastn
          program from R.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## As of Nov 8, 2024, the latest release of V-QUEST is 202430-2:
list_VQUEST_releases()

## Download the V-QUEST 202430-2 germline sequences and store the IG
## sequences for Mouse in a local database:
install_VQUEST_germline_db("202430-2", organism="Mus musculus",
                           db_type="IG", force=TRUE)

## List the local germline databases:
list_germline_dbs()

## Tell igblastn() to use VQUEST-202430-2.Mus_musculus.IG:
use_germline_db("VQUEST-202430-2.Mus_musculus.IG")
}

\keyword{utilities}

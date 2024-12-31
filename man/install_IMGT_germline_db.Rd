\name{install_IMGT_germline_db}

\alias{install_IMGT_germline_db}
\alias{list_IMGT_releases}

\title{Install a germline db from IMGT}

\description{
  The \code{install_IMGT_germline_db()} function downloads IMGT/V-QUEST
  germline sequences (FASTA files) from the IMGT website, and generates
  a small local database with the sequences of the specified organism.
  This local database can then be used later with \code{\link{igblastn}()}.

  Note that \code{install_IMGT_germline_db()} generates the
  local database by performing the instructions provided at
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}.
  In particular, utility tool edit_imgt_file.pl (Perl script included
  in the IgBLAST installation) is invoked to process the sequences.
  This means that Perl must be available on your machine.
}

\usage{
install_IMGT_germline_db(release, organism="Homo sapiens",
                         force=FALSE, ...)

list_IMGT_releases(recache=FALSE)
}

\arguments{
  \item{release}{
    A single string specifying the IMGT/V-QUEST release to get the
    germline sequences from.
    Use \code{list_IMGT_releases()} to list all releases.
  }
  \item{organism}{
    A single string specifying the organism for which to get the germline
    sequences.
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
    \code{list_IMGT_releases()} uses a caching mechanism so that the
    list of IMGT/V-QUEST releases gets downloaded only once from the IMGT
    website during an R session (note that this caching is done in memory
    so it does not persist across sessions).
    Set \code{recache} to \code{TRUE} to force a new download (and
    recaching) of the list of IMGT/V-QUEST releases.
  }
}

\value{
  \code{install_IMGT_germline_db()} returns the name to the local
  germline db as an invisible string.

  \code{list_IMGT_releases()} returns the list of IMGT/V-QUEST releases
  in a character vector. The releases are sorted from newest to oldest
  (latest release is first).
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          program from R. This is the main function in the \pkg{igblastr}
          package.

    \item \code{\link{use_germline_db}} to select the local germline db
          to use with \code{igblastn()}.

    \item The IMGT website: \url{https://www.imgt.org/}.

    \item The IMGT/V-QUEST download site:
          \url{https://www.imgt.org/download/V-QUEST/}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## As of Dec 27, 2024, the latest IMGT/V-QUEST release is 202449-1:
list_IMGT_releases()

## Download germline sequences from IMGT/V-QUEST 202449-1, and store
## the IG sequences for Mouse in a local database:
install_IMGT_germline_db("202449-1", organism="Mus musculus", force=TRUE)

## List the local germline databases:
list_germline_dbs()

## Tell igblastn() to use the newly installed germline db:
use_germline_db("IMGT-202449-1.Mus_musculus.IGH+IGK+IGL")
}

\keyword{utilities}

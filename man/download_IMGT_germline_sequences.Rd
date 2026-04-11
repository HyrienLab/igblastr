\name{download_IMGT_germline_sequences}

\alias{download_IMGT_germline_sequences}
\alias{list_IMGT_releases}
\alias{list_IMGT_organisms}
\alias{IMGT_is_up}

\title{Download germline sequences from IMGT}

\description{
  The \code{download_IMGT_germline_sequences()} function downloads
  germline V/D/J gene allele sequences (nucleotides) from the IMGT/V-QUEST
  download site for a given IMGT/V-QUEST release and organism.
  The sequences are downloaded in several FASTA files, one per
  \emph{IMGT group}.

  CONDITIONS OF USE AND LICENSE: The IMGT data is provided to the academic
  users and NPO's (Not for Profit Organization(s)) under the CC BY-NC-ND 4.0
  license. See \url{https://creativecommons.org/licenses/by-nc-nd/4.0/}.
  Any other use of IMGT material, from the private sector, needs a
  financial arrangement with CNRS.
}

\usage{
download_IMGT_germline_sequences(release, organism="Homo sapiens",
                                 tcr.db=FALSE,
                                 destdir=".", overwrite=FALSE,
                                 recache=FALSE, ...)

## Related utilities:
list_IMGT_releases(recache=FALSE)
list_IMGT_organisms(release)
IMGT_is_up()
}

\arguments{
  \item{release}{
    A single string specifying the IMGT/V-QUEST release to download
    the germline sequences from (or to list the organisms from for
    \code{list_IMGT_organisms()}).
    Use \code{list_IMGT_releases()} to list all releases.
  }
  \item{organism}{
    A single string specifying the latin name of the organism for which
    to download the germline sequences.
  }
  \item{tcr.db}{
    Should the germline sequences to download be B-cell Receptor (a.k.a. BCR)
    sequences or T-cell Receptor (a.k.a. TCR) sequences?

    By default (i.e. when \code{tcr.db} is omitted or set to \code{FALSE}),
    the V/D/J allele sequences for the three BCR loci (IGH, IGK, IGL) will
    be downloaded. This will result in the 7 following FASTA files being
    downloaded:
    \preformatted{        IGHV.fasta, IGHD.fasta, IGHJ.fasta,
        IGKV.fasta, IGKJ.fasta,
        IGLV.fasta, IGLJ.fasta}

    When \code{tcr.db} is set to \code{TRUE}, the V/D/J allele sequences
    for the four TCR loci (TRA, TRB, TRG, TRD) will be downloaded instead.
    This will result in the 10 following FASTA files being downloaded:
    \preformatted{        TRAV.fasta, TRAJ.fasta,
        TRBV.fasta, TRBD.fasta, TRBJ.fasta,
        TRGV.fasta, TRGJ.fasta,
        TRDV.fasta, TRDD.fasta, TRDJ.fasta}

    Note that:
    \itemize{
      \item The names of the downloaded files follow the <IMGT group>.fasta
            convention, where <IMGT group> is one of the 7 \emph{IMGT groups}
	    with the "IG" prefix or one of the 10 \emph{IMGT groups} with
            the "TR" prefix.
      \item IMGT doesn't always provide germline sequences for
            all \emph{IMGT groups}. This means that, depending on the
            IMGT/V-QUEST release and organism, some FASTA files can be missing.
    }
  }
  \item{destdir}{
    A single string that is the path to the directory where the FASTA files
    containing the germline sequences should be downloaded.
  }
  \item{overwrite}{
    \code{download_IMGT_germline_sequences()} won't replace existing
    files, unless \code{overwrite} is set to \code{TRUE}.
  }
  \item{recache}{
    \code{download_IMGT_germline_sequences()} and \code{list_IMGT_releases()}
    both use a caching mechanism for the data that they download from
    IMGT/V-QUEST. In the case of \code{download_IMGT_germline_sequences()},
    the downloaded files are cached on disk so the caching is persistent
    across R sessions. In the case of \code{list_IMGT_releases()}, the
    caching is in memory so does not persist across sessions.

    Set \code{recache} to \code{TRUE} to force a new download (and recaching)
    of the data.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?\link[utils]{download.file}} in the \pkg{utils} package for
    more information.
  }
}

\value{
  \code{download_IMGT_germline_sequences()} returns the names of the
  downloaded FASTA files in an invisible character vector.

  \code{list_IMGT_releases()} returns the list of IMGT/V-QUEST releases
  in a character vector. The releases are sorted from newest to oldest
  (latest release is first).

  \code{list_IMGT_organisms()} returns the list of organisms included
  in the specified IMGT/V-QUEST release in a character vector.

  \code{IMGT_is_up()} returns \code{TRUE} or \code{FALSE}, indicating
  whether the IMGT website at \url{https://www.imgt.org} is up and running
  or down.
}

\seealso{
  \itemize{
    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item \code{\link{install_custom_germline_db}} to install a germline db
          from user-supplied gene allele sequences.

    \item \code{\link{download_OGRDB_germline_sequences}} to download
          germline sequences from OGRDB.

    \item The IMGT website: \url{https://www.imgt.org/}.

    \item The IMGT/V-QUEST download site:
          \url{https://www.imgt.org/download/V-QUEST/}.
  }
}

\examples{
if (IMGT_is_up()) {
  ## As of April 10, 2026, the latest IMGT/V-QUEST release is 202614-2:
  list_IMGT_releases()

  list_IMGT_organisms("202614-2")

  ## Download Mouse BCR germline gene allele sequences from IMGT/V-QUEST
  ## 202614-2 to a temporary directory:
  dir.create(destdir <- tempfile("fasta_dir_"))
  download_IMGT_germline_sequences("202614-2", organism="Mus musculus",
                                   destdir=destdir)

  ## List the downloaded files:
  list.files(destdir)  # 7 FASTA files

  ## Remove temporary directory:
  unlink(destdir, recursive=TRUE)

  ## Download Mouse TCR germline gene allele sequences from IMGT/V-QUEST
  ## 202614-2 to a temporary directory:
  dir.create(destdir <- tempfile("fasta_dir_"))
  download_IMGT_germline_sequences("202614-2", organism="Mus musculus",
                                   tcr.db=TRUE, destdir=destdir)

  ## List the downloaded files:
  list.files(destdir)  # 10 FASTA files

  ## Remove temporary directory:
  unlink(destdir, recursive=TRUE)
}
}

\keyword{utilities}

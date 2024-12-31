\name{igblastn}

\alias{igblastn}
\alias{igblastn_help}

\title{BLAST for Ig and TCR sequences}

\description{
  WORK IN PROGRESS!
  [Please use \url{https://github.com/hpages/igblastr/issues} to report bugs,
  provide feedback, request features, etc...]

  \code{igblastn()} runs the igblastn program from R. This is the main
  function in the \pkg{igblastr} package.
}

\usage{
igblastn(query, outfmt="AIRR", organism="auto", ...,
         show.in.browser=FALSE, show.command.only=FALSE)

igblastn_help(long.help=FALSE, show.in.browser=FALSE)
}

\arguments{
  \item{query}{
    A single string containing the path to a FASTA file with your
    query sequences. Alternatively \code{query} can be a named
    \link[Biostrings]{DNAStringSet} object.
  }
  \item{outfmt}{
    TODO
  }
  \item{organism}{
    \code{"auto"} (the default), or the organism for your query sequences.
    Supported organisms include human, mouse, rat, rabbit and
    rhesus_monkey for Ig and human and mouse for TCR.
    Use \code{\link{list_igblast_organisms}()} to obtain this list
    programmatically.

    Note that by default (i.e. when \code{organism} is omitted or set to
    \code{"auto"}), \code{igblastn()} infers the organism from the name
    of the germline db in use.
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

    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item \code{\link{use_germline_db}} to select the local germline db
          to use with \code{igblastn()}.

    \item \code{\link{use_c_region_db}} to select the local C-region db
          to use with \code{igblastn()}.

    \item \link[Biostrings]{DNAStringSet} objects implemented in the
          \pkg{Biostrings} package.

    \item \code{\link{list_igblast_organisms}()} to list the organisms
          supported by IgBLAST.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## File 'catnap_bnabs.fasta' contains a small BCR sequencing data set
## made of 1000 heavy- and light-chains associated with bnAbs downloaded
## from the CATNAP database:
catnap_bnabs <- system.file(package="igblastr", "extdata",
                            "catnap_bnabs.fasta")

## Install Human germline db from IMGT (Perl required!):
db_name <- install_IMGT_germline_db("202449-1", "Homo_sapiens", force=TRUE)

## Select germline db to use with igblastn():
use_germline_db(db_name)

## Select C-region db to use with igblastn():
use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")

## Run igblastn:
igblastn(catnap_bnabs, show.in.browser=TRUE)
}

\keyword{manip}

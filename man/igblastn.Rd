\name{igblastn}

\alias{print.igblastn_raw_output}
\alias{igblastn}
\alias{igblastn_help}

\title{BLAST for Ig and TCR sequences}

\description{
  The \code{igblastn()} function is a wrapper to the \code{igblastn}
  \emph{standalone executable} included in IgBLAST. This is the main
  function in the \pkg{igblastr} package.
}

\usage{
igblastn(query, outfmt="AIRR",
         germline_db_V="auto", germline_db_V_seqidlist=NULL,
         germline_db_D="auto", germline_db_D_seqidlist=NULL,
         germline_db_J="auto", germline_db_J_seqidlist=NULL,
         organism="auto", c_region_db="auto", auxiliary_data=NULL,
         ...,
         out=NULL, parse.out=TRUE,
         show.in.browser=FALSE, show.command.only=FALSE)

igblastn_help(long.help=FALSE, show.in.browser=FALSE)
}

\arguments{
  \item{query}{
    A single string that is the path to the FASTA file containing
    the query sequences. Alternatively \code{query} can be a named
    \link[Biostrings]{DNAStringSet} object.
  }
  \item{outfmt}{
    One of \code{"AIRR"}, \code{3}, \code{4}, \code{7}, or \code{19}.
    \code{"AIRR"} is the default and is an alias for \code{19}.
    \code{outfmt} can also be a string describing a customized format 7
    e.g. \code{"7 qseqid sseqid pident nident length score"}.
    See \code{?\link{list_outfmt7_specifiers}} for more information about
    customizing format 7.
  }
  \item{germline_db_V}{
    \code{"auto"} (the default), or the path to a V-region db.

    Note that by default (i.e. when \code{germline_db_V} is omitted or set
    to \code{"auto"}), \code{igblastn()} uses the V-region db that belongs
    to the local germline db currently selected.
    See \code{?\link{use_germline_db}} for how to select the local
    germline db to use with \code{igblastn()}.
  }
  \item{germline_db_D}{
    Same as \code{germline_db_V} but for the D-region db.
  }
  \item{germline_db_J}{
    Same as \code{germline_db_V} but for the J-region db.
  }
  \item{germline_db_V_seqidlist,
        germline_db_D_seqidlist,
        germline_db_J_seqidlist}{
    Restrict search of germline database to list of sequence ids.
    Each list of sequence ids can be specified either as a character vector
    of ids or as the path to a file containing the sequence ids (one id per
    line). In the latter case, a \code{file} object must be passed to the
    \code{germline_db_V_seqidlist}, \code{germline_db_D_seqidlist}, or
    \code{germline_db_J_seqidlist} argument. The \code{file} object will
    typically be constructed with something like
    \code{file("path/to/some/file")}.
  }
  \item{organism}{
    \code{"auto"} (the default), or the organism associated with the query
    sequences. Supported organisms include human, mouse, rat, rabbit and
    rhesus_monkey. Use \code{\link{list_igblast_organisms}()} to obtain
    this list programmatically.

    Note that by default (i.e. when \code{organism} is omitted or set
    to \code{"auto"}), \code{igblastn()} infers the organism from the
    name of the local germline db currently selected.
    See \code{?\link{use_germline_db}} for how to select the local
    germline db to use with \code{igblastn()}.
  }
  \item{c_region_db}{
    \code{"auto"} (the default), \code{NULL}, or the path to a C-region db.

    Note that by default (i.e. when \code{c_region_db} is omitted or set
    to \code{"auto"}), \code{igblastn()} uses the local C-region db
    currently selected.
    See \code{?\link{use_c_region_db}} for how to select the local
    C-region db to use with \code{igblastn()}.
  }
  \item{auxiliary_data}{
    \code{NULL} (the default), or the path to a file containing the
    coding frame start positions for the sequences in the J-region db.

    IMPORTANT NOTE:

    Unfortunately, \code{igblastn()} can emit a significant number
    of the following warning when the \code{auxiliary_data} argument
    is not specified:
    \preformatted{    Warning: Auxilary data file could not be found}
    Note that you can use \code{\link{get_igblast_auxiliary_data}()} to
    obtain the path to the auxiliary data included in the IgBLAST
    installation used by the \pkg{igblastr} package for a given organism,
    and pass the returned path to the \code{auxiliary_data} argument.
    This will silence the annoying warnings. However, be aware that this
    data is NOT guaranteed to be compatible with the germline db that you
    are using. See \code{?\link{get_igblast_auxiliary_data}} for the details.

    If you supply auxiliary data that is not compatible with your germline
    dbs, then \code{igblastn()} will NOT return proper frame status or CDR3
    information (other returned information will still be correct).
  }
  \item{...}{
    Extra arguments to be passed to the \code{igblastn}
    \emph{standalone executable}. The list of valid arguments can be
    displayed with \code{igblastn_help()}.

    Note that the argument/value pairs must be passed to the \code{igblastn()}
    function in the usual R fashion. For example, what would be passed
    as \code{-num_alignments 1 -num_threads 8} when invoking the
    \code{igblastn} \emph{standalone executable} in a terminal should
    be passed as \code{num_alignments_V=1, num_threads=8} when
    calling the \code{igblastn()} function:
    \preformatted{    igblastn(query, num_alignments_V=1, num_threads=8)}
  }
  \item{out}{
    \code{NULL} (the default), or the path to the file where the
    \code{igblastn} \emph{standalone executable} should write its output.

    Note that by default (i.e. when \code{out} is omitted or set
    to \code{NULL}), \code{igblastn()} instructs the \code{igblastn}
    \emph{standalone executable} to write its output to a temporary file.
  }
  \item{parse.out}{
    Whether \code{igblastn()} should parse the plain-text output produced
    by the \code{igblastn} \emph{standalone executable} or not, before
    returning it to the user. \code{TRUE} by default.

    If set to \code{FALSE}, then \code{igblastn()} returns the output
    as-is in a character vector, with one line per element in the vector.
    Note that \code{igblastn()} sets the \code{"igblastn_raw_output"} class
    attribute on this character vector, which allows compact display of the
    vector (this is achieved via a dedicated \code{print()} method defined
    in the \pkg{igblastr} package). The class attribute can be dropped with
    \code{unclass()}.
  }
  \item{show.in.browser}{
    For \code{igblastn()}: Whether the output of the \code{igblastn}
    \emph{standalone executable} should also be displayed in the browser
    or not (in addition to being returned by the \code{igblastn()}
    function call). \code{FALSE} by default.

    For \code{igblastn_help()}: Whether the help printed by the
    \code{igblastn} \emph{standalone executable} (when invoked with
    the -h or -help argument) should be displayed in the browser or not.
    \code{FALSE} by default.
  }
  \item{show.command.only}{
    \code{TRUE} or \code{FALSE}. If set to \code{TRUE}, \code{igblastn()}
    won't invoke the \code{igblastn} \emph{standalone executable} and
    instead will display the full command that shows how it would have
    invoked it. Note that the command is also returned in an invisible
    character vector. \code{FALSE} by default.
  }
  \item{long.help}{
    \code{TRUE} or \code{FALSE}. If set to \code{FALSE} (the default),
    the \code{igblastn} \emph{standalone executable} is invoked with
    the -h argument. Otherwise, it's invoked with the -help argument.
  }
}

\value{
  \code{igblastn()} captures the output produced by the \code{igblastn}
  \emph{standalone executable} and returns it as:
  \itemize{
    \item A \link[tibble]{tibble} with 1 row per query sequence if
          \code{outfmt} is \code{"AIRR"} or \code{19} and \code{parse.out}
          is \code{TRUE}.
    \item A nested list with two top-level components (\code{records}
          and \code{footer}) if \code{outfmt} is \code{7} (or a customized
          format 7) and \code{parse.out} is \code{TRUE}.
          See \code{?\link{parse_outfmt7}} for more information.
    \item A character vector with class attribute \code{"igblastn_raw_output"}
          on it in all other cases, that is, if \code{parse.out} is
          \code{FALSE} or \code{outfmt} is \code{3} or \code{4}.
          See the \code{parse.out} argument above for more information.
  }
}

\seealso{
  \itemize{
    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.

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

    \item \code{\link{list_outfmt7_specifiers}} for how to customize
          output format 7.

    \item \code{\link{list_igblast_organisms}()} to list the organisms
          supported by IgBLAST.

    \item \link[Biostrings]{DNAStringSet} objects implemented in the
          \pkg{Biostrings} package.

    \item \link[tibble]{tibble} objects implemented in the
          \pkg{tibble} package.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## ---------------------------------------------------------------------
## Preliminary steps
## ---------------------------------------------------------------------

igblast_info()

## File 'catnap_bnabs.fasta' contains a small BCR sequencing data set
## made of 1000 heavy- and light-chains associated with bnAbs downloaded
## from the CATNAP database:
catnap_bnabs <- system.file(package="igblastr", "extdata",
                            "catnap_bnabs.fasta")

## Install Human germline db from IMGT (Perl required!):
db_name <- install_IMGT_germline_db("202506-1", "Homo_sapiens", force=TRUE)

## Select germline db to use with igblastn():
use_germline_db(db_name)

## Select C-region db to use with igblastn():
use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")

## ---------------------------------------------------------------------
## Call igblastn()
## ---------------------------------------------------------------------

## We don't specify the 'outfmt' argument so output will be in AIRR format.
## Since we're not supplying auxiliary data, igblastn() will emit many
## warnings (see documentation of the 'auxiliary_data' argument above for
## the details). We silence them by wrapping the call to igblastn() in
## suppressWarnings().
out <- suppressWarnings(igblastn(catnap_bnabs))
out

## The result is a tibble with one row per query sequence:
class(out)
dim(out)

## ---------------------------------------------------------------------
## More examples
## ---------------------------------------------------------------------

## See '?parse_outfmt7' for more examples.
}

\keyword{manip}

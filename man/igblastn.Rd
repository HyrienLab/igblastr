\name{igblastn}

\alias{print.igblastn_raw_output}
\alias{igblastn}
\alias{igblastn_help}

\title{BLAST for Ig and TCR sequences}

\description{
  \code{igblastn()} runs the \emph{igblastn} executable from R. This is
  the main function in the \pkg{igblastr} package.
}

\usage{
igblastn(query, outfmt="AIRR", parse.out=TRUE, organism="auto",
         germline_db_V="auto", germline_db_V_seqidlist=NULL,
         germline_db_D="auto", germline_db_D_seqidlist=NULL,
         germline_db_J="auto", germline_db_J_seqidlist=NULL,
         c_region_db="auto", auxiliary_data=NULL, ...,
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
  \item{parse.out}{
    Whether the output of the \emph{igblastn} executable should be parsed
    or not. \code{TRUE} by default. If set to \code{FALSE}, the output is
    a text that \code{igblastn()} returns in a character vector with one
    element per line. Note that the character vector has a class attribute
    \code{"igblastn_raw_output"} on it which allows compact display of the
    vector via a dedicated \code{print()} method. The class attribute can
    be dropped with \code{unclass()}.
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
  \item{c_region_db}{
    \code{"auto"} (the default), \code{NULL}, or the path to a C-region db.

    Note that by default (i.e. when \code{c_region_db} is omitted or set
    to \code{"auto"}), \code{igblastn()} uses the local C-region db
    currently selected.
    See \code{?\link{use_c_region_db}} for how to select the local
    C-region db to use with \code{igblastn()}.
  }
  \item{auxiliary_data}{
    \code{NULL} (the default), or the path to a file containing the coding
    frame start positions for the sequences in the J-region db.

    IMPORTANT NOTE:

    Unfortunately, \code{igblastn()} can emit a significant number
    of the following warning when the \code{auxiliary_data} argument
    is not specified:
    \preformatted{    Warning: Auxilary data file could not be found}
    According to \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html},
    a standard IgBLAST installation -- like the one used by the \pkg{igblastr}
    package -- is expected to include auxiliary data that is specific to a
    particular NCBI or IMGT germline db. Unfortunately this means that this
    data is unlikely to be compatible with the germline db that you are using.

    Note that you can use \code{get_igblast_auxiliary_data()} to obtain
    the path to the auxiliary data included in the IgBLAST installation
    used by the \pkg{igblastr} package for a given organism, and set the
    \code{auxiliary_data} argument to it. This will silence the warnings.
    However, be aware that if this data is not compatible with the
    germline db that you are using, which is very likely to be the case,
    then \code{igblastn()} will NOT return proper frame status or CDR3
    information (other returned information will still be correct).
  }
  \item{...}{
    Extra arguments to be passed to the \emph{igblastn} executable.
  }
  \item{show.in.browser}{
    For \code{igblastn()}: Whether the output of the \emph{igblastn}
    executable should also be displayed in the browser or not (in addition
    to being returned by the \code{igblastn()} function call). \code{FALSE}
    by default.

    For \code{igblastn_help()}: Whether the help printed by the
    \emph{igblastn} executable (when called with the -h or -help argument)
    should be displayed in the browser or not. \code{FALSE} by default.
  }
  \item{show.command.only}{
    \code{TRUE} or \code{FALSE}. If set to \code{TRUE}, \code{igblastn()}
    won't call the \emph{igblastn} executable and instead will display
    the full command that shows how it would have called it. Note that
    the command is also returned in an invisible character vector.
    \code{FALSE} by default.
  }
  \item{long.help}{
    \code{TRUE} or \code{FALSE}. If set to \code{FALSE} (the default),
    the \emph{igblastn} executable is called with the -h argument.
    Otherwise, it's called with the -help argument.
  }
}

\value{
  \code{igblastn()} captures the output of calling the \emph{igblastn}
  executable and returns it as:
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
db_name <- install_IMGT_germline_db("202449-1", "Homo_sapiens", force=TRUE)

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

\name{igblastn}

\alias{print.igblastn_raw_output}
\alias{igblastn}
\alias{igblastn_help}

\title{BLAST for BCR/Ig and TCR sequences}

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
         organism="auto", c_region_db="auto",
         custom_internal_data="auto", auxiliary_data="auto",
         domain_system=c("imgt", "kabat"), ig_seqtype="auto",
         ...,
         out=NULL, parse.out=TRUE,
         show.in.browser=FALSE, show.command.only=FALSE)

igblastn_help(long.help=FALSE, show.in.browser=FALSE)
}

\arguments{
  \item{query}{
    A character vector containing the paths to the input files (FASTA),
    or a \emph{named} \link[Biostrings]{DNAStringSet} object.

    If a character vector, then \code{query} must be of length >= 1 and
    each vector element must be the path to a FASTA file (possibly
    gz-compressed). In the context of IgBLAST, the DNA sequences in the
    FASTA files are referred to as \emph{the query sequences}, and the
    sequence names found in the description lines of the FASTA records
    are referred to as \emph{the query sequence ids}.

    Note that the query sequences are typically (but not always) stored
    in a single file, in which case \code{query} will be a single string.
    If more than one FASTA file is specified via \code{query}, then
    \code{igblastn()} will concatenate all the files together and pass
    the resulting file to the \code{igblastn} \emph{standalone executable}.

    If \code{query} is a \link[Biostrings]{DNAStringSet} object, then it
    must have names on it. These will be considered the query sequence ids.
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

    Note that, by default (i.e. when \code{germline_db_V} is omitted or set
    to \code{"auto"}), \code{igblastn()} uses the V-region db that belongs
    to the cached germline db that is currently selected.

    See \code{?\link{use_germline_db}} for how to select the cached
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
    Restrict search of germline database to list of gene alleles.
    A list of gene alleles can be specified either as a character vector
    of gene allele identifiers (e.g. \code{IGHV3-23*01}, \code{IGHV3-23*04},
    etc...) or as the path to a file containing the identifiers (one
    identifier per line). In the latter case, a \code{file} object must
    be passed to the \code{germline_db_V_seqidlist},
    \code{germline_db_D_seqidlist}, or \code{germline_db_J_seqidlist}
    argument. The \code{file} object will typically be constructed with
    something like \code{file("path/to/some/file")}.
  }
  \item{organism}{
    \code{"auto"} (the default), or the organism associated with the query
    sequences. Supported organisms include human, mouse, rat, rabbit and
    rhesus_monkey. Use \code{\link{list_igblast_organisms}()} to obtain
    this list programmatically.

    Note that, by default (i.e. when \code{organism} is omitted or set
    to \code{"auto"}), \code{igblastn()} infers the organism from the
    name of the cached germline db that is currently selected.

    See \code{?\link{use_germline_db}} for how to select the cached
    germline db to use with \code{igblastn()}.
  }
  \item{c_region_db}{
    \code{"auto"} (the default), \code{NULL}, or the path to a C-region db.

    Note that, by default (i.e. when \code{c_region_db} is omitted or
    set to \code{"auto"}), \code{igblastn()} uses the cached C-region db
    that is currently selected.

    See \code{?\link{use_c_region_db}} for how to select the cached
    C-region db to use with \code{igblastn()}.
  }
  \item{custom_internal_data}{
    \code{"auto"} (the default), or a data.frame containing custom FWR/CDR
    annotation, or the path to a file containing custom FWR/CDR annotation,
    or \code{NULL}.

    By default (i.e. when \code{custom_internal_data} and \code{germline_db_V}
    are both omitted or set to \code{"auto"}), \code{igblastn()} will use:
    \itemize{
        \item the \emph{internal data} included in the cached germline
              db that is currently selected, if it has any;
        \item the organism-specific IgBLAST \emph{internal data} from NCBI
              if the cached germline db that is currently selected does not
              include its own \emph{internal data}.
    }
    Use \code{\link{list_germline_dbs}(with.intdata.only=TRUE)} to list
    the cached germline dbs that include their own \emph{internal data}.

    See \code{?\link{use_germline_db}} for how to select the cached
    germline db to use with \code{igblastn()}.

    Note that when \code{custom_internal_data} is set to \code{NULL},
    \code{igblastn()} always uses the organism-specific IgBLAST
    \emph{internal data} from NCBI.
  }
  \item{auxiliary_data}{
    \code{"auto"} (the default), or a data.frame containing auxiliary data,
    or the path to a file containing auxiliary data (.aux file),
    or \code{NULL}.

    By default (i.e. when \code{auxiliary_data} and \code{germline_db_J}
    are both omitted or set to \code{"auto"}), \code{igblastn()} will use:
    \itemize{
        \item the \emph{auxiliary data} included in the cached germline
              db that is currently selected, if it has any;
        \item the organism-specific IgBLAST \emph{auxiliary data} from NCBI
              if the cached germline db that is currently selected does not
              include its own \emph{auxiliary data}.
    }
    Use \code{\link{list_germline_dbs}(with.auxdata.only=TRUE)} to list
    the cached germline dbs that include their own \emph{auxiliary data}.

    See \code{?\link{use_germline_db}} for how to select the cached
    germline db to use with \code{igblastn()}.

    IMPORTANT NOTES:
    \itemize{
        \item Supplying \emph{auxiliary data} that is not compatible with
              the J gene sequences of the selected germline db can cause
              \code{igblastn()} to return improper frame status or CDR3
              information (other returned information will still be correct).
              See \code{?\link{get_auxdata_path}} for more information.
        \item When \code{auxiliary_data} is set to \code{NULL}, then no
              \emph{auxiliary data} is used. In this case, \code{igblastn()}
              can emit a significant number of the following warning:
              \preformatted{    Warning: Auxilary data file could not be found}
              and various columns of the returned AIRR-formatted
              \link[tibble]{tibble} (e.g. columns \code{vj_in_frame},
              \code{productive}, \code{cdr3}, \code{fwr4}, and others)
              will be filled with \code{NA}s.
      }
  }
  \item{domain_system}{
    Set to \code{"imgt"} or \code{"kabat"}.
  }
  \item{ig_seqtype}{
    Set to \code{"Ig"} or \code{"TCR"} depending on whether the query
    sequences are BCR/Ig or TCR sequences.

    Note that, by default (i.e. when \code{ig_seqtype} is omitted or
    set to \code{"auto"}), the value of \code{ig_seqtype} is inferred
    from the germline loci that appear in the name of the cached
    germline db that is currently selected (this name can be obtained with
    \code{\link{use_germline_db}()}). If these are BCR/Ig germline loci
    (i.e. IGH, IGK, IGL), then the inferred value will be \code{"Ig"}.
    If they are TCR germline loci (TRA, TRB, TRG, TRD), then it will
    be \code{"TCR"}.

    See \code{?\link{use_germline_db}} for how to select the cached
    germline db to use with \code{igblastn()}.
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

    For options that don't require a value (e.g. \code{-extend_align5end},
    \code{-extend_align3end}, \code{-ungapped}, etc...), pass the empty string
    (or a white string) to the argument. For example:
    \preformatted{    igblastn(query, extend_align5end="", extend_align3end="")}
  }
  \item{out}{
    \code{NULL} (the default), or the path to the file where the
    \code{igblastn} \emph{standalone executable} should write its output.

    Note that, by default (i.e. when \code{out} is omitted or set
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
    \emph{standalone executable} should also be displayed in a browser
    or not (in addition to being returned by the \code{igblastn()}
    function call). \code{FALSE} by default.

    For \code{igblastn_help()}: Whether the help printed by the
    \code{igblastn} \emph{standalone executable} (when invoked with
    the -h or -help argument) should be displayed in a browser or not.
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
          See \code{?\link{read_igblastn_fmt7_output}} for more information.
    \item A character vector with class attribute \code{"igblastn_raw_output"}
          on it in all other cases, that is, if \code{parse.out} is
          \code{FALSE} or \code{outfmt} is \code{3} or \code{4}.
          See the \code{parse.out} argument above for more information.
  }
}

\note{
  By default, the NCBI BLAST+ and IgBLAST programs will "call home" to
  report usage when they run on a computer with internet access. See
  \url{https://www.ncbi.nlm.nih.gov/books/NBK569851/} for the details.
  This can induce a significant slowdown in some situations e.g. when
  the \code{igblastn} \emph{standalone executable} is called in a loop
  on a small set of query sequences at each iteration.

  For this reason, the "call home" feature is disabled in \pkg{igblastr}
  by default, unless environment variable BLAST_USAGE_REPORT is set to true.
  See \code{?\link{igblastr_usage_report}} for more information.
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

    \item \code{\link{install_custom_germline_db}} to install a germline db
          from user-supplied gene allele sequences.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item \code{\link{use_c_region_db}} to select the cached C-region db
          to use with \code{igblastn()}.

    \item \code{\link{igbrowser}} to display the annotated sequences
          returned by \code{igblastn()} in a browser.

    \item \code{\link{list_outfmt7_specifiers}} for how to customize
          output format 7.

    \item \code{\link{list_igblast_organisms}} to list the organisms
          supported by IgBLAST.

    \item \code{\link{augment_germline_db}} to add novel gene alleles to
          a germline db.

    \item \link{igblastr_usage_report} to turn "Usage Reporting" on or off.

    \item \link[Biostrings]{DNAStringSet} objects implemented in the
          \pkg{Biostrings} package.

    \item \link[tibble]{tibble} objects implemented in the
          \pkg{tibble} package.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## Access query sequences and select germline and C-region dbs to use
## ---------------------------------------------------------------------

## Files 'heavy_sequences.fasta' and 'light_sequences.fasta' included
## in igblastr contain 250 paired heavy- and light- chain sequences (125
## sequences in each file) downloaded from OAS (the Observed Antibody
## Space database):
filenames <- paste0(c("heavy", "light"), "_sequences.fasta")
query <- system.file(package="igblastr", "extdata", "BCR", filenames)

## Install Human germline db from IMGT:
db_name <- install_IMGT_germline_db("202614-2", "Homo_sapiens",
                                    overwrite=TRUE)

## Select germline db to use with igblastn():
use_germline_db(db_name)

## Select C-region db to use with igblastn():
use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")

## ---------------------------------------------------------------------
## Call igblastn()
## ---------------------------------------------------------------------

## We don't specify the 'outfmt' argument so output will be in AIRR
## format:
AIRR_df <- igblastn(query)
AIRR_df

## The result is a tibble with one row per query sequence:
class(AIRR_df)
dim(AIRR_df)

## You can call igbrowser() on 'AIRR_df' to visualize the annotated
## sequences in a browser. See '?igbrowser'.

## Note that this tibble can easily be converted to an ordinary data.frame
## with 'as.data.frame()', or to a DataFrame with 'as(., "DataFrame")':
as(AIRR_df, "DataFrame")

## To call igblastn() on a subset of the FASTA file, load the file as a
## DNAStringSet object with Biostrings::readDNAStringSet(), then subset
## the object, and finally pass the result of the subsetting operation
## to igblastn():
query_21_30 <- readDNAStringSet(query)[21:30]
query_21_30  # a DNAStringSet object with 10 sequences
igblastn(query_21_30)

## ---------------------------------------------------------------------
## Parallel computing
## ---------------------------------------------------------------------

## The igblastn standalone executable included in IgBLAST supports
## multithreading via command line argument -num_threads. To use it
## in igblastn(), set argument 'num_threads' to the desired value:
AIRR_df2 <- igblastn(query, num_threads=4)
stopifnot(identical(AIRR_df, AIRR_df2))

## Unfortunately, in our experience, using 'num_threads' on Linux
## doesn't achieve any significant speedup.

## Alternatively, one can parallelize execution at the R level using
## standard tools from the parallel or BiocParallel package. In this
## case it's the responsibility of the user to split the input in
## smaller batches and combine the results obtained for each batch:
if (.Platform$OS.type != "windows") {
  library(parallel)
  ## Split 'query' in 10 batches of 25 sequences each:
  batches <- split(readDNAStringSet(query), rep(1:10, each=25))
  AIRR_df3 <- mclapply(seq_along(batches),
                       function(i) igblastn(batches[[i]]),
                       mc.cores=4)
  ## Combine the results:
  AIRR_df3 <- do.call(rbind, AIRR_df3)
  stopifnot(identical(AIRR_df, AIRR_df3))
}

## Note that the actual speedup will depend on many factors like
## hardware, number and size of the batches, number of cores used,
## etc...

## ---------------------------------------------------------------------
## TCR analysis
## ---------------------------------------------------------------------

## NCBI IgBLAST can also be used for TCR sequence analysis, and so does
## igblastn().

## File 'SRR11341217.fasta.gz' included in igblastr contains 10,875 human
## beta chain TCR transcripts running from 5' of reverse transcription
## reaction to beginning of constant region:
filename <- "SRR11341217.fasta.gz"
query <- system.file(package="igblastr", "extdata", "TCR", filename)

## For this example, we're only keeping the first 100 sequences:
query <- head(readDNAStringSet(query), n=100)

## Install Human TCR germline db from IMGT:
db_name <- install_IMGT_germline_db("202614-2", "Homo_sapiens",
                                    tcr.db=TRUE, overwrite=TRUE)

## Select germline db to use with igblastn():
use_germline_db(db_name)

## Select C-region db to use with igblastn():
use_c_region_db("_IMGT.human.TRA+TRB+TRG+TRD.202509")

## Call igblastn(). Note that the 'ig_seqtype' argument will be
## automatically set to "TCR" (see documentation of the 'ig_seqtype'
## argument above in this man page for more information):
AIRR_df <- igblastn(query)

## ---------------------------------------------------------------------
## More examples
## ---------------------------------------------------------------------

## See '?read_igblastn_fmt7_output' for more examples.
}

\keyword{manip}

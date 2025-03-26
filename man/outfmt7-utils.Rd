\name{outfmt7-utils}

\alias{outfmt7-utils}

\alias{list_outfmt7_specifiers}
\alias{print.outfmt7_specifiers}

\alias{qseqid}
\alias{qseqid.query_details}
\alias{summary.query_details}
\alias{print.query_details}
\alias{print.VDJ_rearrangement_summary}
\alias{print.VDJ_junction_details}
\alias{print.subregion_sequence_details}
\alias{print.alignment_summary}
\alias{print.hit_table}
\alias{qseqid.fmt7record}
\alias{print.fmt7record}
\alias{print.fmt7footer}
\alias{parse_outfmt7}

\title{Handle igblastn output format 7}

\description{
  Some utilities to handle igblastn output format 7.
}

\usage{
list_outfmt7_specifiers()

parse_outfmt7(out)
}

\arguments{
  \item{out}{
    The output of \code{igblatsn(query, outfmt=7, parse.out=FALSE, ...)}.
  }
}

\value{
  \code{list_outfmt7_specifiers()} returns the list of format specifiers
  supported by \code{igblastn()} formatting option 7.

  \code{parse_outfmt7(out)} returns the parsed form of \code{out} in a list.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          executable from R. This is the main function in the \pkg{igblastr}
          package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## File 'catnap_bnabs.fasta' contains a small BCR sequencing data set
## made of 1000 heavy- and light-chains associated with bnAbs downloaded
## from the CATNAP database:
catnap_bnabs <- system.file(package="igblastr", "extdata",
                            "catnap_bnabs.fasta")
query <- head(readDNAStringSet(catnap_bnabs), n=10)

db_name <- install_IMGT_germline_db("202506-1", "Homo_sapiens", force=TRUE)
use_germline_db(db_name)
use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")

## ---------------------------------------------------------------------
## FIRST igblastn RUN: GET OUTPUT IN FORMAT 7
## ---------------------------------------------------------------------

## For this first run we specify 'outfmt=7' and 'parse.out=FALSE':
out <- igblastn(query, outfmt=7, parse.out=FALSE)
out  # raw output

out <- parse_outfmt7(out)  # parse the output

## Output contains one record per query sequence:
length(out$records)  # 10

## Each record can have 5 or 6 sections:
##   1. query_details
##   2. VDJ_rearrangement_summary
##   3. VDJ_junction_details
##   4. subregion_sequence_details (can be missing)
##   5. alignment_summary
##   6. hit_table

## Taking a close look at the first record:
rec1 <- out$records[[1]]
rec1

qseqid(rec1)    # query sequence id associated with this record

rec1$hit_table  # data.frame with the standard columns

## ---------------------------------------------------------------------
## SECOND igblastn RUN: GET OUTPUT IN CUSTOMIZED FORMAT 7
## ---------------------------------------------------------------------

## For this second run we request a customized format 7 by supplying
## space delimited format specifiers:
outfmt <- "7 qseqid sseqid pident nident length score"
out <- igblastn(query, outfmt=outfmt)

## Taking a close look at the first record:
rec1 <- out$records[[1]]
rec1$hit_table  # data.frame with the requested columns (+ the
                # automatic "chaintype" column)
}

\keyword{manip}

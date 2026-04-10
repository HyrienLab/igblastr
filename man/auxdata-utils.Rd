\name{auxdata-utils}

\alias{auxdata-utils}
\alias{auxdata_utils}
\alias{auxdata}

\alias{get_auxdata_path}
\alias{load_auxdata}
\alias{translate_J_alleles}
\alias{J_allele_has_stop_codon}
\alias{translate_fwr4}

\alias{get_igblast_auxiliary_data}
\alias{load_igblast_auxiliary_data}

\title{Access, manipulate, and generate IgBLAST auxiliary data}

\description{
  IgBLAST \emph{auxiliary data} is expected to annotate all the known
  germline J gene alleles for a given organism. It is provided by NCBI
  and is typically included in a standard IgBLAST installation.

  The \pkg{igblastr} package provides a small set of utilities to access,
  manipulate, and generate IgBLAST \emph{auxiliary data}.
}

\usage{
## Access auxiliary data:
get_auxdata_path(organism, which=c("live", "original"))
load_auxdata(organism, which=c("live", "original"))

## Manipulate auxiliary data:
translate_J_alleles(J_alleles, auxdata)
J_allele_has_stop_codon(J_alleles, auxdata)
translate_fwr4(J_alleles, auxdata, max.codons=NA)
}

\arguments{
  \item{organism}{
    A single string containing the name of an organism as
    returned by \code{\link{list_igblast_organisms}()}.

    Alternatively, this can be the name of a cached germline db.
    Note that this works only for germline dbs that include their
    own \emph{auxiliary data}. You can use:
    \preformatted{    list_germline_dbs(with.auxdata.only=TRUE)}
    to list them.
    See \code{?\link{list_germline_dbs}} for more information.
  }
  \item{which}{
    By default, \code{get_auxdata_path()} and \code{load_auxdata()}
    access the "live IgBLAST data", that is, the IgBLAST data that
    the user has possibly updated with \code{update_live_igdata()}.
    Depending on whether updates were applied or not, the "live IgBLAST data"
    might differ from the original IgBLAST data.

    Set \code{which} to \code{"original"} if you want to access
    the original IgBLAST data instead.

    See \code{?\link{update_live_igdata}} for more information about
    "live" and "original" IgBLAST data.
  }
  \item{J_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline J
    gene allele sequences.
  }
  \item{auxdata}{
    A data.frame as returned by \code{load_auxdata()} or
    \code{\link{compute_auxdata}()}.
  }
  \item{max.codons}{
    The maximum number of FWR4 codons to translate. By default (i.e. when
    \code{max.codons} is \code{NA}) all the FWR4 codons are translated.
  }
}

\details{
  IgBLAST \emph{auxiliary data} is typically included in a standard IgBLAST
  installation. It's located in the \code{optional_file/} directory which
  is itself a subdirectory of IgBLAST \emph{root directory}.

  The data consists of one tabulated file per organism. Each file indicates
  the germline J gene coding frame start position, the J gene type, and the
  CDR3 end position for all known germline J allele sequences. See
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html} for
  additional details.
}

\value{
  \code{get_auxdata_path()} returns a single string containing
  the path to the \emph{auxiliary data} included in the IgBLAST installation
  used by \pkg{igblastr}, for the specified organism. Not necessarily suitable
  to use with \code{\link{igblastn}()} (see WARNING below).

  \code{load_auxdata()} returns the \emph{auxiliary data} in a data.frame
  with 1 row per supplied J allele sequence and the same columns as the
  data.frame returned by \code{read_auxdata()}.
  See \code{?\link{read_auxdata}} for more information.

  \code{translate_J_alleles()} returns a named character vector with 1
  amino acid sequence per supplied allele. The vector contains an \code{NA}
  for any allele that is not annotated in \code{auxdata} or for which
  \code{auxdata$coding_frame_start} has an \code{NA}. The names on it are
  the names of the supplied alleles.

  \code{J_allele_has_stop_codon()} returns a named logical vector with 1
  value per supplied allele. The vector contains an \code{NA} for any
  allele that is not annotated in \code{auxdata} or for which
  \code{auxdata$coding_frame_start} has an \code{NA}. The names on it are
  the names of the supplied alleles.

  \code{translate_fwr4()} returns a named character vector with 1 amino
  acid sequence per supplied allele. The vector contains an \code{NA}
  for any allele that is not annotated in \code{auxdata} or for which
  \code{auxdata$cdr3_end} has an \code{NA}.
}

\section{WARNING}{
  According to \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
  the \emph{auxiliary data} included in IgBLAST is specific to a particular
  NCBI or IMGT germline db. Unfortunately this means that this data is
  NOT guaranteed to be compatible with the germline db that you will
  use with \code{\link{igblastn}()}. See documentation of the
  \code{auxiliary_data} argument in \code{?\link{igblastn}} for
  more information about this.
}

\seealso{
  \itemize{
    \item \code{\link{compute_auxdata}} to annotate a set of germline J
          gene allele sequences.

    \item \link{intdata_utils} to access and manipulate IgBLAST
          \emph{internal data}.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about IgBLAST \emph{auxiliary data}.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link{translate_codons}} function on which
          \code{translate_J_alleles()} and \code{translate_fwr4()} are based.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## 1. Access and load IgBLAST auxiliary data for a given organism
## ---------------------------------------------------------------------

list_igblast_organisms()

get_auxdata_path("human")
human_auxdata <- load_auxdata("human")
head(human_auxdata)

## Get the auxiliary data included in germline db
## _OGRDB.rhesus_monkey.IGH+IGK+IGL.202602 (this data
## was computed with compute_auxdata()):
db_name <- "_OGRDB.rhesus_monkey.IGH+IGK+IGL.202602"
rhesus_monkey_auxdata <- load_intdata(db_name)
head(rhesus_monkey_auxdata)

## ---------------------------------------------------------------------
## 2. A close look at IgBLAST auxiliary data for rabbit
## ---------------------------------------------------------------------

get_auxdata_path("rabbit")
rabbit_auxdata <- load_auxdata("rabbit")

## It turns out that IgBLAST auxiliary data for rabbit matches exactly
## the set of rabbit germline J alleles available at IMGT:
db_name <- install_IMGT_germline_db("202531-1", "Oryctolagus cuniculus",
                                    without.auxdata=TRUE,
                                    overwrite=TRUE)
J_alleles <- load_germline_db(db_name, region_types="J")
J_alleles  # DNAStringSet object
stopifnot(setequal(names(J_alleles), rabbit_auxdata$allele_name))

## Note that this might change with future IMGT releases.

## Let's put the allele sequences in 'J_alleles' in the same order as
## in 'rabbit_auxdata':
J_alleles <- J_alleles[rabbit_auxdata$allele_name]
stopifnot(identical(names(J_alleles), rabbit_auxdata$allele_name))

## The 'coding_frame_start' column in 'rabbit_auxdata' contains integer
## values that are >= 0 and <= 2. They indicate how many nucleotides
## precede the first codon on each allele sequence. In other words,
## this is the number of nucleotides that we need to trim on the 5'
## end of the germline J allele sequence before we start translating
## it. translate_J_alleles() uses this information to translate the
## DNA sequences in 'J_alleles':
J_aa <- translate_J_alleles(J_alleles, rabbit_auxdata)
J_aa

## No sequence in 'J_aa' should contain the letter "*" which is used
## by translate_J_alleles() to represent a stop codon. However, one
## J allele in IMGT-202531-1.Oryctolagus_cuniculus.IGH+IGK+IGL seems
## to disobey:
has_stop_codon <- grepl("*", J_aa, fixed=TRUE)
rabbit_auxdata[has_stop_codon, ]  # coding_frame_start = 0 for IGKJ1-2*04
J_alleles[has_stop_codon]         # first codon (TGA) is a stop codon
J_aa[has_stop_codon]              # indeed!

## ---------------------------------------------------------------------
## 3. About the "WGXG" and "FGXG" motifs
## ---------------------------------------------------------------------

## The FWR4 region is expected to start with the following amino acid
## motifs (X represents any amino acid):
##   - "WGXG" on the heavy chain
##   - "FGXG" on the light chain

## Let's use translate_fwr4() to extract and translate the first 4
## codons of the FWR4 region:
fwr4_head <- translate_fwr4(J_alleles, rabbit_auxdata, max.codons=4)

## We expect to see the "WGXG" and "FGXG" motifs here, and most of the
## time we do:
has_motif <- grepl("[FW]G.G", fwr4_head)
table(has_motif)

## However, there are a few exceptions:
fwr4_head[!has_motif]

## ---------------------------------------------------------------------
## 4. Compute auxiliary data for a set of J allele sequences
## ---------------------------------------------------------------------

## See '?compute_auxdata' for the details of what compute_auxdata() does.
computed_auxdata <- compute_auxdata(J_alleles)
head(computed_auxdata)

## Alleles for which the CDR3 end could not found:
cdr3_end_not_found <- is.na(computed_auxdata$cdr3_end)
stopifnot(identical(cdr3_end_not_found, !has_motif))
J_alleles[cdr3_end_not_found]
fwr4_head[cdr3_end_not_found]  # déjà vu

## 'computed_auxdata' is in agreement with 'rabbit_auxdata', except for
## the 8 alleles for which the CDR3 end could not be found:
keep_idx <- which(!cdr3_end_not_found)
stopifnot(identical(computed_auxdata[keep_idx, ],
                    rabbit_auxdata[keep_idx, ]))
}

\keyword{utilities}

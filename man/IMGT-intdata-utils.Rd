\name{IMGT-intdata-utils}

\alias{IMGT-intdata-utils}
\alias{IMGT_intdata_utils}
\alias{get_fwrcdr_widths_for_imgt_organism}

\title{Get the FWR/CDR widths for a given IMGT organism}

\description{
  \code{get_fwrcdr_widths_for_imgt_organism()} returns the "real"
  FWR/CDR widths for a given IMGT organism.

  Note that it is used internally by \code{\link{install_IMGT_germline_db}()}
  to compute the \emph{internal data} that gets included in the germline db.
}

\usage{
get_fwrcdr_widths_for_imgt_organism(organism, loci_prefix=c("IG", "TR"))
}

\arguments{
  \item{organism}{
    A single string specifying the latin name of the organism for which
    to get the FWR/CDR widths.
  }
  \item{loci_prefix}{
    Either \code{"IG"} to get the FWR/CDR widths for \emph{IMGT groups}
    IG[HKL]V, or \code{"TR"} to get them for \emph{IMGT groups} TR[ABGD]V.
  }
}

\details{
  For some IMGT organisms, the gapped V allele sequences in the
  IG[HKL]V.fasta or TR[ABGD]V.fasta files provided by IMGT can have
  additional insertions in their FWR/CDR regions w.r.t. to the gaps
  already inserted according to the IMGT unique numbering.

  The additional insertions can be seen by going to the IMGT "Protein
  displays" section at \url{https://www.imgt.org/IMGTrepertoire/Proteins/}
  and clicking on any organism for any \emph{IMGT group} (i.e. any IG[HKL]V
  or TR[ABGD]V group).

  Because of these additional insertions, the widths of the FWR/CDR regions
  (i.e. the number of letters+dots used to represent these regions in
  the gapped amino acid sequences of the V alleles) can be more than the
  default. See \code{?\link{IMGT_DEFAULT_FWRCDR_WIDTHS}} for the default
  widths.

  \code{get_fwrcdr_widths_for_imgt_organism()} returns the "real"
  FWR/CDR widths for a given IMGT organism, that is, the FWR/CDR
  widths that account for the additional insertions. The returned
  vector or matrix is meant to be passed to the \code{fwrcdr_widths}
  argument of \code{\link{compute_V_gene_delineations}()}.
}

\value{
  Can return any of the following:
  \itemize{
    \item An \code{NA} if no IG[HKL]V.fasta (or TR[ABGD]V.fasta) files
          are provided for the specified organism, or the FWR/CDR widths
          are not known.
    \item \code{IMGT_DEFAULT_FWRCDR_WIDTHS} if the IG[HKL]V.fasta (or
          TR[ABGD]V.fasta) files provided by IMGT for the specified
          organism have no additional insertions.
    \item A named numeric vector similar to \code{IMGT_DEFAULT_FWRCDR_WIDTHS}
          but with possibly some of the values greater than the values
          in \code{IMGT_DEFAULT_FWRCDR_WIDTHS} if the IG[HKL]V.fasta (or
          TR[ABGD]V.fasta) files provided by IMGT for the specified
          organism have additional insertions.
    \item A 5-column numeric matrix with 1 row per \emph{IMGT group}.
          The colnames on the matrix are fwr1, cdr1, fwr2, cdr2, and fwr3.
          The rownames on the matrix are the \emph{IMGT groups} for which
          the corresponding IG[HKL]V.fasta (or TR[ABGD]V.fasta) file exists
          and the FWR/CDR widths are known.
  }
}

\seealso{
  \itemize{
    \item \url{https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html}
          for the role of gaps in the germline V gene allele sequences
          provided by IMGT.

    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.

    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.
  }
}

\examples{
## No additional insertions in the IG[HKL]V.fasta or TR[ABGD]V.fasta
## files provided by IMGT for human:
get_fwrcdr_widths_for_imgt_organism("Homo sapiens")
get_fwrcdr_widths_for_imgt_organism("Homo sapiens", loci_prefix="TR")

## For rhesus monkey:
## - IGHV.fasta has 2 additional insertions in the FWR1
## - IGKV.fasta has 1 additional insertion in the FWR1
## - IGLV.fasta has 1 additional insertion in the FWR1 and 2 in the FWR2
get_fwrcdr_widths_for_imgt_organism("Macaca_mulatta")

## For rat:
## - IGHV.fasta has 1 insertion in the CDR1 and 1 in the FWR3
## - IGKV.fasta and IGLV.fasta have no additional insertions
get_fwrcdr_widths_for_imgt_organism("Rattus_norvegicus")

## IMGT provides no TR[ABGD]V.fasta files for rat:
get_fwrcdr_widths_for_imgt_organism("Rattus_norvegicus", loci_prefix="TR")

## For zebrafish (TR loci):
## - TRAV.fasta has 1 additional insertion in the FWR1, 1 in the FWR2
##   (between pos 46 & 47), 2 in the CDR2, and 2 in the FWR3 (between
##   pos 74 & 75 and between pos 84 & 85)
## - TRBV.fasta, TRGV.fasta, and TRDV.fasta are not provided
get_fwrcdr_widths_for_imgt_organism("Danio_rerio", loci_prefix="TR")

## For sheep (TR loci):
## - TRAV.fasta has 3 additional insertions in the FWR2
## - TRBV.fasta has no additional insertions
## - TRGV.fasta is not provided
## - TRDV.fasta has 1 additional insertion in the CDR1
get_fwrcdr_widths_for_imgt_organism("Ovis_aries", loci_prefix="TR")

## Sanity checks:
stopifnot(
  identical(
    get_fwrcdr_widths_for_imgt_organism("Homo sapiens"),
    IMGT_DEFAULT_FWRCDR_WIDTHS
  ),
  identical(
    get_fwrcdr_widths_for_imgt_organism("Homo sapiens", loci_prefix="TR"),
    IMGT_DEFAULT_FWRCDR_WIDTHS
  ),
  identical(
    get_fwrcdr_widths_for_imgt_organism("Rattus_norvegicus"),
    rbind(IGHV=c(fwr1=26L, cdr1=13L, fwr2=17L, cdr2=10L, fwr3=40L),
          IGKV=IMGT_DEFAULT_FWRCDR_WIDTHS,
          IGLV=IMGT_DEFAULT_FWRCDR_WIDTHS)
  ),
  identical(
    get_fwrcdr_widths_for_imgt_organism("Rattus_norvegicus",
                                        loci_prefix="TR"),
    NA
  ),
  identical(
    get_fwrcdr_widths_for_imgt_organism("Danio_rerio", loci_prefix="TR"),
    rbind(TRAV=c(fwr1=27L, cdr1=12L, fwr2=18L, cdr2=12L, fwr3=41L))
  ),
  identical(
    get_fwrcdr_widths_for_imgt_organism("Ovis_aries", loci_prefix="TR"),
    rbind(TRAV=c(fwr1=26L, cdr1=12L, fwr2=20L, cdr2=10L, fwr3=39L),
          TRBV=IMGT_DEFAULT_FWRCDR_WIDTHS,
          TRDV=c(fwr1=26L, cdr1=13L, fwr2=17L, cdr2=10L, fwr3=39L))
  )
)
}

\keyword{utilities}

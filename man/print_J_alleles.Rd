\name{print_J_alleles}

\alias{print_J_alleles}

\title{Fancy display of J alleles sequences}

\description{
  Use \code{print_J_alleles()} to display a set of J allele sequences
  where the sequences are colored and justified w.r.t. their CDR3/FWR4
  junction.
}

\usage{
print_J_alleles(J_alleles, auxdata, translate=FALSE,
                igblast_organism=NA, filler=".", cdr3fwr4_sep=" ")
}

\arguments{
  \item{J_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline J
    gene allele sequences.

    Alternatively, this can be the name of a cached germline db, in
    which case the \code{auxdata} argument doesn't need to (and
    shouldn't) be supplied.
    Note that this works only for germline dbs that include their
    own \emph{auxiliary data}. You can use:
    \preformatted{    list_germline_dbs(with.auxdata.only=TRUE)}
    to list them.
    See \code{?\link{list_germline_dbs}} for more information.
  }
  \item{auxdata}{
    A data.frame as returned by \code{\link{load_auxdata}()},
    \code{\link{compute_auxdata}()}, or
    \code{\link{compute_germline_db_auxdata}()}, that contains
    annotations for the alleles in \code{J_alleles}.
  }
  \item{translate}{
    By default, \code{print_J_alleles()} displays the nucleotide sequences.
    Set \code{translate} to \code{TRUE} to display the amino acid sequences
    instead.
  }
  \item{igblast_organism}{
    If the alleles in \code{J_alleles} are from an "IgBLAST organism" (see
    \code{?\link{list_igblast_organisms}}), then the name of the organism
    can be passed to \code{igblast_organism}.
    In this case \code{print_J_alleles()} will display whether each allele
    is also annotated in the \emph{auxiliary data} provided by IgBLAST for
    this organism.
  }
  \item{filler}{
    The character to use to pad the sequences on both ends. This padding
    achieves the horizontal alignment of all sequences w.r.t. the
    CDR3/FWR4 junction.
  }
  \item{cdr3fwr4_sep}{
    The string to use to separate the CDR3 portions of the sequences from
    their FWR4 portions.
  }
}

\value{
  \code{print_J_alleles()} returns the CDR3 and FWR4 portions of the
  J allele sequences in an invisible \link[S4Vectors]{DataFrame} with
  1 row per J allele.
}

\seealso{
  \itemize{
    \item \link{J_alleles_inspect} for basic inspection of J allele sequences.

    \item \code{\link{load_auxdata}} to access IgBLAST \emph{auxiliary data}.

    \item \code{\link{print_gapped_V_alleles}} to display a set of gapped
          V allele sequences where the sequences are colored and split by
          FWR/CDR regions.

    \item \link{V_alleles_inspect} for basic inspection of V allele sequences.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \link[Biostrings]{DNAStringSet} and \link[Biostrings]{AAStringSet}
          objects in the \pkg{Biostrings} package.

    \item The \code{\link{translate_codons}} function which is used internally
          by \code{print_J_alleles()} when \code{translate=TRUE}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## 1. A QUICK LOOK AT THE GERMLINE J GENE ALLELES PROVIDED BY OGRDB
## ---------------------------------------------------------------------

list_germline_dbs()

## --- Human germline J gene alleles ---

db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"

## The 'also_in_IgBLAST_auxdata' column indicates whether an allele is
## also annotated in the "auxiliary data" provided by IgBLAST for human
## (i.e. in the optional_file/human_gl.aux file shipped with IgBLAST):
print_J_alleles(db_name)
print_J_alleles(db_name, translate=TRUE)

## Note that the above shows that:
## - All FWR4 regions from the IGH locus (heavy chain) start with
##   the "WGXG" motif.
## - All FWR4 regions from the IGK and IGL loci (light chain) start
##   with the "FGXG" motif.

## --- Mouse germline J gene alleles ---

db_name <- "_OGRDB.mouse.MSM_MsJ.IGH+IGK+IGL.202603"

## None of the mouse J alleles from OGRDB are annotated in the "auxiliary
## data" provided by IgBLAST for mouse (i.e. in optional_file/mouse_gl.aux).
print_J_alleles(db_name, translate=TRUE)

## The "auxiliary data" included in this germline db is provided by OGRDB.
## In fact, the "auxiliary data" included in **all** the preinstalled OGRDB
## dbs is provided by OGRDB.

## ---------------------------------------------------------------------
## 2. A QUICK LOOK AT THE GERMLINE J GENE ALLELES PROVIDED BY IMGT
## ---------------------------------------------------------------------

## --- For "IgBLAST organisms" ---

## Mouse:
db_name <- install_IMGT_germline_db("202631-1", "Mus_musculus",
                                    overwrite=TRUE)
print_J_alleles(db_name, translate=TRUE)

## Note that:
## - The annotations for the 3 alleles not annotated in
##   IgBLAST "auxiliary data" for mouse were computed by
##   install_IMGT_germline_db().
## - The CDR3 of allele IGHJ3*02 contains a stop codon (indicated by
##   character "*").

## Rat:
db_name <- suppressWarnings(
  install_IMGT_germline_db("202631-1", "Rattus_norvegicus", overwrite=TRUE)
)
print_J_alleles(db_name, translate=TRUE)

## Notes:
## - The ?> and <? tags around a sequence indicate an allele for which
##   the CDR3/FWR4 junction is unknown.
## - Even though allele IGKJ3*01 is annotated in the "auxiliary data"
##   provided by IgBLAST for rat (in optional_file/rat_gl.aux), its
##   cdr3_end value is missing:
igblast_auxdata <- load_auxdata("rat")
igblast_auxdata
subset(igblast_auxdata, is.na(cdr3_end))

## --- For non-IgBLAST organisms ---

## Gorilla:
db_name <- install_IMGT_germline_db("202631-1", "Gorilla_gorilla_gorilla",
                                    overwrite=TRUE)
print_J_alleles(db_name, translate=TRUE)

## Note that, in this case:
## - There's no 'also_in_IgBLAST_auxdata' column.
## - install_IMGT_germline_db() took care of computing the annotations
##   for all the alleles.

## Pig:
db_name <- install_IMGT_germline_db("202631-1", "Sus_scrofa",
                                    overwrite=TRUE)
print_J_alleles(db_name, translate=TRUE)

## ---------------------------------------------------------------------
## 3. MORE ABOUT print_J_alleles()
## ---------------------------------------------------------------------

## print_J_alleles() also accepts an arbitrary set of J allele sequences
## and their annotations (via arguments 'J_alleles' and 'auxdata',
## respectively):

db_name <- "IMGT-202631-1.Mus_musculus.IGH+IGK+IGL"
J_alleles <- load_germline_sequences(db_name, region_types="J")
auxdata <- load_auxdata(db_name)
print_J_alleles(J_alleles, auxdata, igblast_organism="mouse")

## Note that the above is equivalent to 'print_J_alleles(db_name)'.

## print_J_alleles() returns a DataFrame with the CDR3 and FWR4
## portions of the J allele sequences:

DF <- suppressMessages(print_J_alleles(db_name, translate=TRUE))
DF

## ---------------------------------------------------------------------
## 4. FWR4 CONSENSUS SEQUENCES
## ---------------------------------------------------------------------

## The above can be used to determine the consensus amino acid sequence
## for mouse FWR4 regions:

fwr4_aa <- setNames(DF$fwr4, DF$allele_name)
fwr4_aa <- fwr4_aa[width(fwr4_aa) != 0L]
IGHJ_fwr4_aa <- fwr4_aa[substr(names(fwr4_aa), 1, 3) == "IGH"]
IGKJ_fwr4_aa <- fwr4_aa[substr(names(fwr4_aa), 1, 3) == "IGK"]
IGLJ_fwr4_aa <- fwr4_aa[substr(names(fwr4_aa), 1, 3) == "IGL"]
consensusString(IGHJ_fwr4_aa)  # WGQGT?VTVSS  (starts with "WGXG" motif)
consensusString(IGKJ_fwr4_aa)  # FG?GTKLEIK   (starts with "FGXG" motif)
consensusString(IGLJ_fwr4_aa)  # FGGGT?LTVL   (starts with "FGXG" motif)

## See '?consensusString' in the Biostrings package for more information
## about the consensusString() function.

## Let's take a look at the FWR4 sequences for human TCR J alleles:
db_name <- suppressWarnings(
  install_IMGT_germline_db("202631-1", "Homo_sapiens", tcr.db=TRUE,
                           overwrite=TRUE)
)
print_J_alleles(db_name, translate=TRUE)

DF <- suppressMessages(print_J_alleles(db_name, translate=TRUE))
fwr4_aa <- setNames(DF$fwr4, DF$allele_name)
fwr4_aa <- fwr4_aa[width(fwr4_aa) != 0L]
TRAJ_fwr4_aa <- fwr4_aa[substr(names(fwr4_aa), 1, 3) == "TRA"]
TRBJ_fwr4_aa <- fwr4_aa[substr(names(fwr4_aa), 1, 3) == "TRB"]
TRGJ_fwr4_aa <- fwr4_aa[substr(names(fwr4_aa), 1, 3) == "TRG"]
TRDJ_fwr4_aa <- fwr4_aa[substr(names(fwr4_aa), 1, 3) == "TRD"]
consensusString(TRAJ_fwr4_aa)  # FG?GT?L?V?P
consensusString(TRBJ_fwr4_aa)  # FG?GTRLTVL
consensusString(TRGJ_fwr4_aa)  # FG?GT?LIVTSP
consensusString(TRDJ_fwr4_aa)  # FGKGT?L?VEP

## Note that, unlike with BCR J alleles where the FWR4 is expected to
## start with the "WGXG" and "FGXG" motifs on the heavy- and light-chain,
## respectively, for human TCR J alleles it's expected to start with
## the "FGXG" motif regardless of the locus.
}

\keyword{utilities}

\name{install_custom_germline_db}

\alias{install_custom_germline_db}

\title{Install a germline db from user-supplied gene allele sequences}

\description{
  The \code{install_custom_germline_db()} function creates a germline
  database from a set of user-supplied FASTA files containing germline
  V/D/J gene allele sequences (nucleotides), and installs it in
  \pkg{igblastr}'s persistent cache. The new germline database can
  then be used later with \code{\link{igblastn}()}.
}

\usage{
install_custom_germline_db(db_name, fasta_dir, tcr.db=FALSE,
                           loci="auto", gapped=FALSE, with.intdata=FALSE,
                           with.auxdata=FALSE, imgt.fasta=FALSE,
                           disambiguate.allele.names=FALSE,
                           overwrite=FALSE, verbose=FALSE)
}

\arguments{
  \item{db_name}{
    A single string specifying the name of the germline db to install.
    The string \emph{must} start with "cus" and cannot contain whitespace
    characters.
  }
  \item{fasta_dir}{
    The path to the directory containing the FASTA files to import in
    the database. Note that the files to import \emph{must} be named:
    \itemize{
      \item IGH[VDJ].fasta, IGK[VJ].fasta, and IGL[VJ].fasta for
            a BCR germline loci database;
      \item TRA[VJ].fasta, TRB[VDJ].fasta, TRG[VJ].fasta, and
            TRD[VDJ].fasta for a TCR germline loci database.
    }
    See the \code{tcr.db} argument below for more details.

    Note that the FASTA header lines are expected to contain the allele
    names:
    \itemize{
      \item either on their own if the header lines have no \code{|},
      \item or in the 2nd field if the header lines have fields
            separated with \code{|}.
    }
  }
  \item{tcr.db}{
    Should the database to install be populated with allele sequences
    from the BCR (B-cell Receptor) or TCR (T-cell Receptor) germline loci?

    The BCR germline loci are: IGH, IGK, IGL.

    The TCR germline loci are: TRA, TRB, TRG, TRD.

    By default, the V/D/J allele sequences from the BCR germline loci
    will be used. More precisely, the 7 following FASTA files will be
    expected to be present in \code{fasta_dir}:
    \preformatted{        IGHV.fasta, IGHD.fasta, IGHJ.fasta,
        IGKV.fasta, IGKJ.fasta,
        IGLV.fasta, IGLJ.fasta}
    The files that are effectively present will be imported in the database,
    with a warning if some of them are missing.

    If \code{tcr.db} is set to \code{TRUE}, then the V/D/J allele sequences
    from the TCR germline loci will be used instead. More precisely, the 10
    following FASTA files will be expected to be present in \code{fasta_dir}:
    \preformatted{        TRAV.fasta, TRAJ.fasta,
        TRBV.fasta, TRBD.fasta, TRBJ.fasta,
        TRGV.fasta, TRGJ.fasta,
        TRDV.fasta, TRDD.fasta, TRDJ.fasta}
    The files that are effectively present will be imported in the database,
    with a warning if some of them are missing.
  }
  \item{loci}{
    By default, the database to install will be populated with the allele
    sequences from all the BCR or TCR loci.
    However, if you want to restrict the database to specific loci,
    you can use the \code{loci} argument to specify these loci.
    The subset of loci can be specified as a character vector with one
    element per locus (e.g. \code{"IGH"} or \code{c("TRA", "TRB")}), or
    as a \code{+}-separated list in a single string (e.g. \code{"TRA+TRB"}).
  }
  \item{gapped}{
    Set to \code{TRUE} if the supplied V allele sequences are \emph{gapped},
    in which case the gaps will be removed when the sequences get imported
    in the database.
  }
  \item{with.intdata}{
    If the supplied V allele sequences are \emph{gapped}, you can
    set \code{with.intdata} to \code{TRUE} to have the \emph{internal data}
    (i.e. the FWR/CDR boundaries on the V alleles) computed and added to
    the new db.

    Note that the presence of \emph{internal data} in the new db will be
    indicated in the \code{intdata} column of the data.frame returned by
    \code{\link{list_germline_dbs}()}.

    If a cached germline db includes its own \emph{internal data}, then
    \code{\link{igblastn}()} will use it (as \emph{custom internal data})
    instead of the \emph{internal data} shipped with IgBLAST.
    See documentation of the \code{custom_internal_data} argument in
    \code{?\link{igblastn}} for more information.

    Also see \code{?\link{intdata_utils}} for more information about IgBLAST
    \emph{internal data}.
  }
  \item{with.auxdata}{
    COMING SOON...
  }
  \item{imgt.fasta}{
    COMING SOON...
  }
  \item{disambiguate.allele.names}{
    Note that \emph{repeated} alleles (i.e. alleles with identical
    ungapped DNA sequences \emph{and} names) are dropped. More precisely,
    only the first allele in each group of \emph{repeated} alleles is kept.

    If, after dropping the \emph{repeated} alleles, the names of the
    remaining alleles are not unique, then an error will be raised, unless
    \code{disambiguate.allele.names} is set to \code{TRUE}, in which case
    the \emph{ambiguous} allele names will be disambiguated by adding a
    suffix to them.

    Finally, note that we only look at \emph{repeated} alleles or
    \emph{ambiguous} allele names within the same region type e.g.
    within the union of IGHV.fasta, IGKV.fasta, and IGLV.fasta (V regions),
    or within the union of TRBD.fasta and TRDD.fasta (D regions).
    But we never look at \emph{repeated} alleles or \emph{ambiguous}
    allele names across IGHV.fasta and IGHJ.fasta, or across TRBV.fasta
    and TRBD.fasta (these are very very unlikely to occur anyways).
  }
  \item{overwrite}{
    Set to \code{TRUE} if a germline db with the specified name is
    already installed in \pkg{igblastr}'s persistent cache, in which
    case the existing db will be replaced with the new one.
  }
  \item{verbose}{
    Set to \code{TRUE} to have the function display some details about
    its internal operations.
  }
}

\value{
  \code{install_custom_germline_db()} returns the name to the newly
  installed germline db as an invisible string.
}

\note{
  \code{install_custom_germline_db()} creates the local
  database by performing the instructions provided at
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}.
}

\seealso{
  \itemize{
    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{load_germline_db}} to load the nucleotide sequences
          of the gene regions stored in a cached germline db.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item \code{\link{load_intdata}} to access IgBLAST \emph{internal data}.

    \item \code{\link{download_IMGT_germline_sequences}} to download
          germline sequences from IMGT.

    \item \code{\link{download_OGRDB_germline_sequences}} to download
          germline sequences from OGRDB.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## ---------------------------------------------------------------------
## A. WITH UNGAPPED GERMLINE V GENE ALLELE SEQUENCES
## ---------------------------------------------------------------------

## In this admittedly artificial example, we're going to create and
## install a germline db from the gene allele sequences in built-in
## germline db _OGRDB.human.IGH+IGK+IGL.202410. More precisely, we're
## going to:
##   1. Dump built-in germline db _OGRDB.human.IGH+IGK+IGL.202410
##      in 7 FASTA files.
##   2. Use these FASTA files to create and install a new germline db.
##   3. Check that the new db contains the same sequences as the
##      original db.

## --- 1. Dump _OGRDB.human.IGH+IGK+IGL.202410 in 7 FASTA files ---

dump_germline_db <- function(db_name, destdir) {
  for (region_type in c("V", "D", "J")) {
    alleles <- load_germline_db(db_name, region_types=region_type)
    for (locus in c("IGH", "IGK", "IGL")) {
      loc_alleles <- alleles[grep(locus, names(alleles))]
      if (length(loc_alleles) == 0L)
        next
      fasta_file <- paste0(locus, region_type, ".fasta")
      writeXStringSet(loc_alleles, file.path(destdir, fasta_file))
    }
  }
}

db_name0 <- "_OGRDB.human.IGH+IGK+IGL.202410"
dir.create(exA_fasta_dir <- tempfile())

dump_germline_db(db_name0, exA_fasta_dir)

list.files(exA_fasta_dir)  # 7 FASTA files

## --- 2. Use FASTA files to install new germline db ---

exA_db_name <- "cus.exA.human.IGH+IGK+IGL"
install_custom_germline_db(exA_db_name, exA_fasta_dir)
list_germline_dbs()

## --- 3. Check that new and original dbs have the same sequences ---

stopifnot(identical(
    as.character(load_germline_db(db_name0)),
    as.character(load_germline_db(exA_db_name))
))

## ---------------------------------------------------------------------
## B. WITH GAPPED GERMLINE V GENE ALLELE SEQUENCES
## ---------------------------------------------------------------------

## In this example, we're going to download human germline gene allele
## **gapped** sequences from AIRR-community/OGRDB and use them to create
## and install a new germline db. More precisely, we're going to:
##   1. Use download_OGRDB_germline_sequences() to download the
##      gapped sequences for the following "Reference Sets" from
##      AIRR-community/OGRDB:
##        - IGH_VDJ version 9
##        - IGKappa_VJ version 4
##        - IGLambda_VJ version 3
##      Each "Reference Set" corresponds to a BCR germline locus. For
##      each set, download_OGRDB_germline_sequences() will download the
##      FASTA file provided by AIRR-community/OGRDB that contains the
##      gapped allele sequences for that locus, and split the file into
##      one FASTA file per region type.
##   2. Use these FASTA files to create and install a new germline db.
##   3. Check that the new db contains the same sequences
##      as built-in germline db _OGRDB.human.IGH+IGK+IGL.202410.

## --- 1. Download human gapped sequences from AIRR-community/OGRDB ---

dir.create(exB_fasta_dir <- tempfile())
germline_sets <- c(IGH_VDJ=9, IGKappa_VJ=4, IGLambda_VJ=3)

download_OGRDB_germline_sequences("Homo sapiens", germline_sets,
                                  destdir=exB_fasta_dir)

list.files(exB_fasta_dir)  # 7 FASTA files

## Note that the V allele sequences have gaps:
IGKV_alleles <- readDNAStringSet(file.path(exB_fasta_dir, "IGKV.fasta"))
as.character(IGKV_alleles[1:3])

## --- 2. Use FASTA files to install new germline db ---

exB_db_name <- "cus.exB.human.IGH+IGK+IGL"
install_custom_germline_db(exB_db_name, exB_fasta_dir, gapped=TRUE)
list_germline_dbs()

## --- 3. Check that the new db contains the same sequences ---
## ---          as _OGRDB.human.IGH+IGK+IGL.202410          ---

## The Reference Sets we used to populate the new db (IGH_VDJ version 9,
## IGKappa_VJ version 4, IGLambda_VJ version 3) are the same that were
## used to populate built-in germline db _OGRDB.human.IGH+IGK+IGL.202410.
## Let's check that the two dbs contain the same sequences:
stopifnot(identical(
    as.character(load_germline_db(db_name0)),
    as.character(load_germline_db(exB_db_name))
))

## ---------------------------------------------------------------------
## C. CREATE AND INSTALL GERMLINE DB WITH COMPUTED "INTERNAL DATA"
## ---------------------------------------------------------------------

## The gaps that germline gene allele sequence providers like IMGT and
## AIRR-community/OGRDB inject in the V allele sequences reflect their
## FWR/CDR boundaries. When the 'with.intdata' argument is set to TRUE,
## install_custom_germline_db() takes advantage of this to compute and
## add the "internal data" associated with the V alleles to the new db:
exC_db_name <- "cus.exC.human.IGH+IGK+IGL"
install_custom_germline_db(exC_db_name, exB_fasta_dir, gapped=TRUE,
                           with.intdata=TRUE)
list_germline_dbs()  # 'intdata' col indicates presence of "internal data"

## Note that the "internal data" in a germline db can be loaded with
## load_intdata() (see '?load_intdata' for more information):
intdataC <- load_intdata(exC_db_name)
head(intdataC)

## Check that the "internal data" in the new db is the same as in
## _OGRDB.human.IGH+IGK+IGL.202410:
stopifnot(identical(intdataC, load_intdata(db_name0)))

## ---------------------------------------------------------------------
## Remove the 3 custom dbs
## ---------------------------------------------------------------------

rm_germline_db(exA_db_name)
rm_germline_db(exB_db_name)
rm_germline_db(exC_db_name)
}

\keyword{utilities}

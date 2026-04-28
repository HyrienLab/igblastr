\name{augment_germline_db}

\alias{augment_germline_db}
\alias{augment_germline_db_V}
\alias{augment_germline_db_D}
\alias{augment_germline_db_J}

\title{Add novel gene alleles to a germline db}

\description{
  WARNING: Some shortcomings were identified in the design of the
  \code{augment_germline_db_[VDJ]()} functions so they are now deprecated.
  Please do not use them. They will be replaced with a better alternative
  in future versions of the package.

  Three functions to add novel V, D, or J gene alleles to a germline db.

  Note that these functions can also be used to combine germline databases
  from two different organisms. See "COMBINE GERMLINE DATABASES FROM TWO
  ORGANISMS" in the Examples section below for how to do this.
}

\usage{
augment_germline_db_V(db_name, novel_alleles,
                      destdir=".", overwrite=FALSE, verbose=FALSE)

augment_germline_db_D(db_name, novel_alleles,
                      destdir=".", overwrite=FALSE, verbose=FALSE)

augment_germline_db_J(db_name, novel_alleles,
                      destdir=".", overwrite=FALSE, verbose=FALSE)
}

\arguments{
  \item{db_name}{
    A single string that is the name of the cached germline db that
    contains the set of gene alleles to augment.
    Use \code{\link{list_germline_dbs}()} to list the cached germline dbs.

    The exact function used (i.e. \code{augment_germline_db_V()},
    \code{augment_germline_db_D()}, or \code{augment_germline_db_J()})
    determines the set of alleles to augment (i.e. alleles from the V, D,
    or J region).
  }
  \item{novel_alleles}{
    A single string that is the path to a FASTA file (possibly gz-compressed)
    where the novel alleles are stored.

    Alternatively, the novel alleles can be supplied as a \emph{named}
    \link[Biostrings]{DNAStringSet} object.
  }
  \item{destdir}{
    A single string that is the path to the "destination directory", that is,
    the directory where the augmented V-, D-, or J-region db is to be created.
    This directory will be created if it doesn't exist already. Note that, by
    default, the augmented region db will be created in the current directory.
  }
  \item{overwrite}{
    If the "destination directory" already contains a V-, D-, or J-region db,
    should it be overwritten?
  }
  \item{verbose}{
    Set to \code{TRUE} to have the function display some details about its
    internal operations.
  }
}

\value{
  These functions don't return anything (invisible \code{NULL}).
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
\dontrun{
if (!has_igblast()) install_igblast()

query <- system.file(package="igblastr", "extdata",
                     "BCR", "heavy_sequences.fasta")

use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")

## ---------------------------------------------------------------------
## USE HUMAN GERMLINE DATABASE FROM AIRR
## ---------------------------------------------------------------------

use_germline_db("_OGRDB.human.IGH+IGK+IGL.202410")

AIRR_df <- igblastn(query)

## ---------------------------------------------------------------------
## ADD NOVEL V ALLELES
## ---------------------------------------------------------------------

## 'fake_human_V_alleles.fasta' contains made-up novel V alleles:
## - 2 novel alleles for gene IGHV1-8: IGHV1-8*fake1, IGHV1-8*fake2
## - 1 novel allele for gene IGHV4-61: IGHV4-61*fake
my_novel_V_alleles <- system.file(package="igblastr", "extdata",
                                  "novel_germline_alleles",
                                  "fake_human_V_alleles.fasta")

## Take a quick look at these novel V alleles:
readDNAStringSet(my_novel_V_alleles)

## Create a new V germline database that combines the V alleles
## from _OGRDB.human.IGH+IGK+IGL.202410 with our novel V alleles:
my_Vdb_path <- file.path(tempdir(), "myVdb")
augment_germline_db_V("_OGRDB.human.IGH+IGK+IGL.202410",
                      my_novel_V_alleles,
                      destdir=my_Vdb_path)

## To use this new augmented V germline database with igblastn(),
## supply its path via the 'germline_db_V' argument:
AIRR_df2 <- igblastn(query, germline_db_V=my_Vdb_path)

## ---------------------------------------------------------------------
## A QUICK COMPARISON BETWEEN 'AIRR_df' AND 'AIRR_df2'
## ---------------------------------------------------------------------

## Index of rows where "v_call" has changed between 'AIRR_df'
## and 'AIRR_df2':
idx <- which(AIRR_df$v_call != AIRR_df2$v_call)
idx  # 2 rows

AIRR_df[idx, c("v_call", "v_cigar", "v_identity")]

AIRR_df2[idx, c("v_call", "v_cigar", "v_identity")]

## Besides these 2 rows, all the other rows are the same:
stopifnot(all.equal(AIRR_df[-idx, ], AIRR_df2[-idx, ]))

## ---------------------------------------------------------------------
## COMBINE GERMLINE DATABASES FROM TWO ORGANISMS
## ---------------------------------------------------------------------

## The augment_germline_db_[VDJ]() functions can be used to combine
## germline databases from two different organisms. This can be useful
## for example when working with BCR sequences from mice that have been
## engineered to have both mouse and some human immunoglobulin genes.
##
## To create a hybrid human/mouse V germline database, we can either:
##
## (1) Add all (or a subset of) mouse V alleles to all human V alleles.
##     This is done by extracting mouse V germline allele sequences from
##     a cached germline database and using them to augment a cached
##     germline database for human.
##
## (2) Add all (or a subset of) human V alleles to all mouse V alleles.
##     This is done by extracting human V germline allele sequences from
##     a cached germline database and using them to augment a cached
##     germline database for mouse.
##
## Note that:
## - We can choose to subset or not the V germline allele sequences
##   extracted from one V germline database before adding them to the
##   other V germline database.
## - The two approaches above are equivalent if we don't subset, that
##   is, if we combine **all** human V alleles with **all** mouse V
##   alleles.
## - However if our engineered mice only have a small known subset of
##   human immunoglobulin genes (e.g. IGHV1-2), then we might want to
##   create a hybrid human/mouse germline database that only adds the
##   human alleles for genes IGHV1-2 to the mouse V alleles. In this
##   case we need to use (2).

## Let's do (2):

db_name1 <- "_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410"
db_name2 <- "_OGRDB.human.IGH+IGK+IGL.202410"

## Extract human V germline alleles:
human_V_alleles <- load_germline_db(db_name2, "V")

## Subset to keep only alleles for genes IGHV1-2:
idx <- grep("^IGHV[12]", names(human_V_alleles))
human_V12_alleles <- human_V_alleles[idx]

## Create a new V germline database that combines the mouse V
## alleles from 'db_name1' with the alleles in 'human_V12_alleles':
engmouseVdb_path <- file.path(tempdir(), "engmouseVdb")
augment_germline_db_V(db_name1, human_V12_alleles,
                      destdir=engmouseVdb_path)
}

## Then, assuming that 'query' contains BCR sequences from the
## engineered mice:
\dontrun{
  use_germline_db(db_name1)
  use_c_region_db("_IMGT.mouse.IGH.202509")
  igblastn(query, germline_db_V=engmouseVdb_path, ...)
}

## Note that, by default, the mouse-only D and J databases that we
## selected above with 'use_germline_db(db_name1)' are being used.
## If we also want to create hybrid D and J databases, we need
## to repeat the above steps for each of them. Then we need to
## specify the paths to the 3 hybrid databases when we call igblastn():
\dontrun{
  igblastn(query, germline_db_V=engmouseVdb_path,
                  germline_db_D=engmouseDdb_path,
                  germline_db_J=engmouseJdb_path,
                  ...)
}
}

\keyword{utilities}

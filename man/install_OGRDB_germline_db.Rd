\name{install_OGRDB_germline_db}

\alias{install_OGRDB_germline_db}
\alias{make_intdata_files_from_ogrdb_jsons}
\alias{make_auxdata_files_from_ogrdb_jsons}

\title{Install a germline db from OGRDB germline sets}

\description{
  The \code{install_OGRDB_germline_db()} function creates and installs
  an IgBLAST-compatible germline database from germline V, D, and J gene
  sequences available at OGRDB:
  \url{https://ogrdb.airr-community.org/}

  The new database gets installed in \pkg{igblastr}'s persistent cache and
  can be used later with \code{\link{igblastn}()}.

  If you use OGRDB in your research, please cite:
  \preformatted{  The current landscape of adaptive immune receptor genomic
  and repertoire data: OGRDB and VDJbase
  Lees, Peres, Klein, Amos et al., Nucleic Acids Research, November 2025.
  https://doi.org/10.1093/nar/gkaf1094}
}

\usage{
install_OGRDB_germline_db(organism, germline_sets,
                          suffix="", source_set=FALSE,
                          overwrite=FALSE, verbose=FALSE, ...)
}

\arguments{
  \item{organism}{
    A single string specifying the latin name of the organism for which
    to install the germline db.
    See \url{https://ogrdb.airr-community.org/germline_sets} for the
    organisms that are currently available at OGRDB.
  }
  \item{germline_sets, source_set}{
    See \code{?\link{download_OGRDB_germline_sequences}} for the
    documentation of these arguments.
  }
  \item{suffix}{
    The suffix to add to the name of the installed database.
  }
  \item{overwrite}{
    Set to \code{TRUE} to reinstall if the requested database is already
    installed.
  }
  \item{verbose}{
    Set to \code{TRUE} to have the function display some details about
    its internal operations.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?\link[utils]{download.file}} in the \pkg{utils} package for
    more information.
  }
}

\details{
  \code{install_OGRDB_germline_db()} installs the newly created blast
  database in \pkg{igblastr}'s persistent cache. The database can then
  be used later with \code{\link{igblastn}()}.

  The following naming scheme is used to form the name of the installed
  database:
  \preformatted{    OGRDB.<organism>.<loci>.<suffix>}
  where:
  \enumerate{
    \item <organism> is the latin name (a.k.a. binomial name) of the
          organism with all spaces replaced with underscores (_).
          For example Homo_sapiens or Macaca_mulatta.
    \item <loci> is a string obtained by concatenating the germline
          loci together separated with the + sign. For example IGH+IGK+IGL.
          Note that, for some organisms, OGRDB supports only a subset
          of the 3 IG loci. For example, for Oncorhynchus mykiss, only
          IGH is supported. This will be automatically reflected in the
          name of the installed germline db.
    \item <suffix> is the string supplied via the \code{suffix} argument.
  }
}

\value{
  \code{install_OGRDB_germline_db()} returns the name to the newly
  installed germline db as an invisible string.
}

\seealso{
  \itemize{
    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item \code{\link{install_custom_germline_db}} to install a germline db
          from user-supplied gene allele sequences.

    \item \code{\link{download_OGRDB_germline_sequences}} to download
          germline sequences from OGRDB.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item The OGRDB website: \url{https://ogrdb.airr-community.org/}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## See https://ogrdb.airr-community.org/germline_sets for which organisms
## are available at OGRDB.

## ---------------------------------------------------------------------
## Install White-footed mouse germline db from OGRDB
## ---------------------------------------------------------------------

## In this first example, we show how to install an IgBLAST-compatible
## germline database for Peromyscus leucopus (White-footed mouse). See
## https://ogrdb.airr-community.org/germline_sets/Peromyscus%20leucopus
## for which germline sets are available for Peromyscus leucopus.
germline_sets <- c(IGH_VDJ=2, IGK_VJ=2, IGL_VJ=1)
db_name <- install_OGRDB_germline_db("Peromyscus leucopus", germline_sets)

## List the cached germline databases:
list_germline_dbs()

## Call use_germline_db() on the db name to use it with igblastn():
use_germline_db(db_name)

list_germline_dbs()

rm_germline_db(db_name)

## ---------------------------------------------------------------------
## Install BALB/c/ByJ mouse germline db from OGRDB
## ---------------------------------------------------------------------

## In this second example, we show how to install an IgBLAST-compatible
## germline database for mouse substrain BALB/c/ByJ.
## See https://ogrdb.airr-community.org/germline_sets/Mus%20musculus for
## which germline sets are available for mouse.
germline_sets <- c(`BALB/c IGH`=3,
                   `BALB/c/ByJ IGKV`=1, `BALB/c/ByJ IGLV`=2,
                   `IGKJ (all strains)`=1, `IGLJ (all strains)`=1)
db_name <- install_OGRDB_germline_db("Mus musculus", germline_sets)

## List the cached germline databases:
list_germline_dbs()

## Call use_germline_db() on the db name to use it with igblastn():
use_germline_db(db_name)

list_germline_dbs()

rm_germline_db(db_name)
}

\keyword{utilities}

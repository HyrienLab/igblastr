\name{list_germline_dbs}

\alias{list_germline_dbs}
\alias{print.germline_dbs_df}
\alias{rm_germline_db}

\title{List cached germline dbs}

\description{
  List the \emph{cached germline dbs}, that is, the germline databases
  currently installed in \pkg{igblastr}'s persistent cache.
}

\usage{
## List the cached germline dbs:
list_germline_dbs(builtin.only=FALSE,
                  with.intdata.only=FALSE,
                  with.auxdata.only=FALSE,
                  names.only=FALSE, long.listing=FALSE)

## Remove a germline db from igblastr's persistent cache:
rm_germline_db(db_name)
}

\arguments{
  \item{builtin.only}{
    By default \code{list_germline_dbs()} returns the list of all cached
    germline dbs, including built-in germline dbs. Set \code{builtin.only}
    to \code{TRUE} to return only the list of built-in germline dbs.
    Note that built-in dbs are prefixed with an underscore (\code{_}).
  }
  \item{with.intdata.only}{
    Set \code{with.intdata.only} to \code{TRUE} to return only the list
    of germline dbs that include their own \emph{internal data} (i.e. the
    FWR/CDR boundaries on the V alleles). Note that whether a germline
    db includes its own \emph{internal data} or not is reported in the
    \code{intdata} column.

    If a cached germline db includes its own \emph{internal data}, then
    \code{\link{igblastn}()} will use it (as \emph{custom internal data})
    instead of the \emph{internal data} shipped with IgBLAST.
    See documentation of the \code{custom_internal_data} argument in
    \code{?\link{igblastn}} for more information.
  }
  \item{with.auxdata.only}{
    Set \code{with.auxdata.only} to \code{TRUE} to return only the list
    of germline dbs that include their own \emph{auxiliary data} (i.e.
    annotations of the J alleles). Note that whether a germline
    db includes its own \emph{auxiliary data} or not is reported in the
    \code{auxdata} column.

    If a cached germline db includes its own \emph{auxiliary data}, then
    \code{\link{igblastn}()} will use it instead of the \emph{auxiliary data}
    shipped with IgBLAST.
    See documentation of the \code{auxiliary_data} argument in
    \code{?\link{igblastn}} for more information.
  }
  \item{names.only}{
    By default \code{list_germline_dbs()} returns the list of cached
    germline dbs in a data.frame with one db per row. Set \code{names.only}
    to \code{TRUE} to return the db names only, in which case they are
    returned in a character vector.
  }
  \item{long.listing}{
    \code{TRUE} or \code{FALSE}. If set to \code{TRUE}, then
    \code{list_germline_dbs()} returns a named list with one
    list element per germline db. Each list element is an integer
    matrix that indicates the number of germline gene allele sequences
    per locus and region type.

    Ignored if \code{names.only} is set to \code{TRUE}.
  }
  \item{db_name}{
    A single string specifying the name of the germline db to remove
    from the cache. This cannot be a built-in db.
  }
}

\details{
  \subsection{Cached germline dbs and C-region dbs}{
    The \pkg{igblastr} package provides utility functions to manage the
    \emph{cached germline dbs} and \emph{cached C-region dbs} to use with
    \code{\link{igblastn}()}.

    Terminology used across \pkg{igblastr} documentation:
    \itemize{
      \item A \emph{cached germline db} contains the nucleotide sequences of
            the germline V, D, and J gene alleles for a given organism.
      \item A \emph{cached C-region db} contains the nucleotide sequences of
            the constant gene regions for a given organism.
    }

    This man page documents utilities that manage cached germline dbs.

    See \code{?\link{list_c_region_dbs}} for utilities that manage
    cached C-region dbs.
  }
  \subsection{Built-in dbs}{
    The \pkg{igblastr} package comes with preinstalled cached germline and
    C-region dbs that are called \emph{built-in dbs}.

    Built-in dbs have their name prefixed with an underscore (\code{_}).

    Note that additional dbs can easily be installed with functions
    like \code{\link{install_IMGT_germline_db}} or
    \code{\link{install_custom_germline_db}}.
  }
  \subsection{AIRR-community/OGRDB built-in cached dbs}{
    Note that the built-in germline dbs starting with \code{_AIRR}
    are made of the AIRR-community/OGRDB germline sets available at
    \url{https://ogrdb.airr-community.org/germline_sets/Homo\%20sapiens}
    for human and at
    \url{https://ogrdb.airr-community.org/germline_sets/Mus\%20musculus}.
    for mouse.

    Each AIRR db is populated with the latest germline datasets that were
    available at AIRR-community/OGRDB at the time indicated by the date (in
    YYYYMM format) embedded in the db name.

    The AIRR dbs with the \code{.src} suffix contain the \emph{Source Sets}.
    The AIRR dbs without the \code{.src} suffix contain the \emph{Reference
    Sets}.
    See \url{https://github.com/HyrienLab/igblastr/tree/devel/inst/extdata/germline_sets/AIRR/human/202410/README.md}
    for more information.

    The AIRR-community/OGRDB maintainers recommend to use the
    \emph{Reference Sets} for AIRR-seq analysis.
    See for example \url{https://ogrdb.airr-community.org/germline_set/75}
  }
}

\value{
  \code{list_germline_dbs()} returns the list of all cached germline dbs
  in a data.frame with one db per row (if \code{names.only} is \code{FALSE},
  which is the default), or in a character vector (if \code{names.only}
  is \code{TRUE}). Columns V, D, J in the data.frame indicate the number
  of germline gene allele sequences for each region in each db.

  \code{rm_germline_db()} does not return anything (invisible \code{NULL}).
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item \code{\link{install_custom_germline_db}} to install a germline db
          from user-supplied gene allele sequences.

    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

    \item \code{\link{reset_germline_dbs}} to reset the cached germline dbs.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## Get list of built-in germline dbs only.
list_germline_dbs(builtin.only=TRUE)
list_germline_dbs(builtin.only=TRUE, names.only=TRUE)  # db names only

## Get list of germline dbs that include their own internal data:
list_germline_dbs(with.intdata.only=TRUE)
list_germline_dbs(with.intdata.only=TRUE, names.only=TRUE)

## Long listing:
list_germline_dbs(long.listing=TRUE)
list_germline_dbs(with.intdata.only=TRUE, long.listing=TRUE)

if (IMGT_is_up()) {
  ## Install Mouse germline db from IMGT:
  install_IMGT_germline_db("202614-2", "Homo sapiens", overwrite=TRUE)

  list_germline_dbs()  # all germline dbs

  ## Select germline db to use with igblastn():
  db_name <- "IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL"
  use_germline_db(db_name)  # select germline db to use

  use_germline_db()  # get current selection
}
}

\keyword{utilities}

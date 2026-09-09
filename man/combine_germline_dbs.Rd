\name{combine_germline_dbs}

\alias{combine_germline_dbs}
\alias{combine_c_region_dbs}

\title{Combine two existing germline dbs}

\description{
  Use the \code{combine_germline_dbs()} function to combine two existing
  germline databases, typically (but not necessarily) from two different
  organisms, into a new one.
  The new germline database gets added to \pkg{igblastr}'s persistent
  cache and can then be used later with \code{\link{igblastn}()}.

  Similarly, the \code{combine_c_region_dbs()} function is provided to
  combine two existing C-region databases.
}

\usage{
combine_germline_dbs(db_name, db_name1, db_name2, suffix1, suffix2,
                     overwrite=FALSE, verbose=FALSE)

combine_c_region_dbs(db_name, db_name1, db_name2, suffix1, suffix2,
                     overwrite=FALSE, verbose=FALSE)
}

\arguments{
  \item{db_name}{
    A single string specifying the name of the combined database.
    The string \emph{must} start with "comb" and cannot contain whitespace
    characters.
  }
  \item{db_name1, db_name2}{
    COMING SOON...
  }
  \item{suffix1, suffix2}{
    COMING SOON...
  }
  \item{overwrite}{
    Set to \code{TRUE} if a database with the specified name is already
    installed in \pkg{igblastr}'s persistent cache, in which
    case the existing database will be replaced with the new one.
  }
  \item{verbose}{
    Set to \code{TRUE} to have the function display some details about
    its internal operations.
  }
}

\value{
  \code{combine_germline_dbs()} and \code{combine_c_region_dbs()} both
  return the name to the newly installed database as an invisible string.
}

\seealso{
  \itemize{
    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item \code{\link{list_c_region_dbs}} to list the cached C-region dbs.

    \item \code{\link{use_c_region_db}} to select the cached C-region db
          to use with \code{igblastn()}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## TO BE EXPANDED SOON...

db_name1 <- install_IMGT_germline_db("202631-1", "Homo_sapiens",
                                     overwrite=TRUE)
db_name2 <- install_IMGT_germline_db("202631-1", "Mus_musculus",
                                     overwrite=TRUE)
db_name <- "comb.IMGT-202631-1.human+mouse.IGH+IGK+IGL"
combine_germline_dbs(db_name, db_name1, db_name2,
                     suffix1="_Hs", suffix2="_Mm", overwrite=TRUE)

rm_germline_db(db_name)
}

\keyword{utilities}

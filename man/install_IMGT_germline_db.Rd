\name{install_IMGT_germline_db}

\alias{install_IMGT_germline_db}
\alias{install_IMGT_c_region_db}

\title{Install a germline db from IMGT}

\description{
  The \code{install_IMGT_germline_db()} function downloads germline
  V/D/J gene allele sequences from the IMGT/V-QUEST download site for
  a given IMGT/V-QUEST release and organism, and stores them in a
  local IgBLAST-compatible germline database. This local database gets
  installed in \pkg{igblastr}'s persistent cache. It can then be used
  later with \code{\link{igblastn}()}.

  CONDITIONS OF USE AND LICENSE: The IMGT data is provided to the academic
  users and NPO's (Not for Profit Organization(s)) under the CC BY-NC-ND 4.0
  license. See \url{https://creativecommons.org/licenses/by-nc-nd/4.0/}.
  Any other use of IMGT material, from the private sector, needs a
  financial arrangement with CNRS.
}

\usage{
install_IMGT_germline_db(release, organism="Homo sapiens", tcr.db=FALSE,
                         loci="auto",
                         without.intdata=FALSE, without.auxdata=FALSE,
                         overwrite=FALSE, verbose=FALSE, ...)

## Can only be used to "subset" a built-in C-region databse (see
## advanced example in Examples section below):
install_IMGT_c_region_db(organism, loci,
                         disambiguate.allele.names=FALSE,
                         overwrite=FALSE, verbose=FALSE)
}

\arguments{
  \item{release, organism, tcr.db}{
    See \code{?\link{download_IMGT_germline_sequences}} for the
    documentation of these arguments.
    Note that, in addition to the organism latin name (e.g. "Homo sapiens"),
    \code{install_IMGT_c_region_db()} also accepts the common name of the
    organism (e.g. "human").
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
  \item{without.intdata}{
    By default, \code{install_IMGT_germline_db()} will also compute and
    add the \emph{internal data} (i.e. the FWR/CDR boundaries on the V
    alleles) to the new db. Note that \code{install_IMGT_germline_db()}
    uses the \emph{gapped} V allele sequences provided by IMGT for that.
    You can skip this step by setting \code{without.intdata} to \code{TRUE}.

    The presence of \emph{internal data} in the new db will be indicated
    in the \code{intdata} column of the data.frame returned by
    \code{\link{list_germline_dbs}()}.

    If a cached germline db includes its own \emph{internal data}, then
    \code{\link{igblastn}()} will use it (as \emph{custom internal data})
    instead of the \emph{internal data} shipped with IgBLAST.
    See documentation of the \code{custom_internal_data} argument in
    \code{?\link{igblastn}} for more information.

    Also see \code{?\link{intdata_utils}} for more information about IgBLAST
    \emph{internal data}.
  }
  \item{without.auxdata}{
    COMING SOON...
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
  \item{disambiguate.allele.names}{
    See documentation of the \code{disambiguate.allele.names}
    in \code{?\link{install_custom_germline_db}}.
  }
}

\details{
  The following naming scheme is used to form the name of the installed
  database:
  \preformatted{    IMGT-<release>.<organism>.<loci>}
  where:
  \enumerate{
    \item <release> is the IMGT/V-QUEST release e.g. 202518-3 or 202449-1.
          Use \code{list_IMGT_releases()} to get the list of releases
          currently available at IMGT/V-QUEST.
    \item <organism> is the latin name (a.k.a. binomial name) of the
          organism with all spaces replaced with underscores (_).
          For example Homo_sapiens or Macaca_mulatta.
          Use \code{list_IMGT_organisms("<release>")} to get the list
          of organisms included in a given IMGT/V-QUEST release.
          Note that, starting with release 202405-2, IMGT/V-QUEST
          provides BCR and TCR germline gene allele sequences for mouse
          strain C57BL6J (Mus_musculus_C57BL6J).
    \item <loci> is a string obtained by concatenating the germline
          loci together separated with the + sign. For example IGH+IGK+IGL
          or TRA+TRB+TRG+TRD.
          The list of loci depends on whether the germline gene allele
          sequences for BCR or TCR were requested. See \code{tcr.db} argument
          above for more information.
          Note that for some IMGT/V-QUEST releases/organisms, only a subset
          of the loci are available. For example, in release 202343-3,
          the only TCR germline loci available for Mus_musculus_C57BL6J
          are TRA and TRB. This will be automatically reflected in the
          name of the installed germline db.
  }
}

\value{
  \code{install_IMGT_germline_db()} returns the name to the newly
  installed germline db as an invisible string.

  \code{install_IMGT_c_region_db()} returns the name to the newly
  installed C-region db as an invisible string.
}

\note{
  \code{install_IMGT_germline_db()} creates the local
  database by performing the instructions provided at
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}.
}

\seealso{
  \itemize{
    \item \code{\link{install_custom_germline_db}} to install a germline db
          from user-supplied gene allele sequences.

    \item \code{\link{list_IMGT_releases}} to list IMGT/V-QUEST releases.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item The IMGT website: \url{https://www.imgt.org/}.

    \item The IMGT/V-QUEST download site:
          \url{https://www.imgt.org/download/V-QUEST/}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

if (IMGT_is_up()) {
  ## -------------------------------------------------------------------
  ## BASIC EXAMPLES
  ## -------------------------------------------------------------------

  ## As of March 26, 2025, the latest IMGT/V-QUEST release is 202518-3:
  list_IMGT_releases()

  list_IMGT_organisms("202518-3")

  ## Download human BCR germline gene allele sequences from IMGT/V-QUEST
  ## release 202518-3, and store them in a cached germline database:
  install_IMGT_germline_db("202518-3", organism="Homo sapiens",
                           overwrite=TRUE)

  ## List the cached germline databases:
  list_germline_dbs()

  ## Select newly installed germline db to use with igblastn():
  use_germline_db("IMGT-202518-3.Homo_sapiens.IGH+IGK+IGL")

  ## Download human TCR germline gene allele sequences from IMGT/V-QUEST
  ## release 202518-3, and store them in a cached germline database:
  install_IMGT_germline_db("202518-3", organism="Homo sapiens",
                           tcr.db=TRUE, overwrite=TRUE)
  list_germline_dbs()

  ## -------------------------------------------------------------------
  ## ADVANCED EXAMPLES
  ## -------------------------------------------------------------------

  ## Install an IMGT database for a subset of TCR loci:
  install_IMGT_germline_db("202518-3", organism="Homo sapiens",
                           loci="TRA+TRB", overwrite=TRUE)
  list_germline_dbs()

  ## Install the corresponding C-region database:
  install_IMGT_c_region_db("human", "TRA+TRB", overwrite=TRUE)
  list_c_region_dbs()

  ## Note that install_IMGT_c_region_db() can only be used to "subset"
  ## a built-in C-region databse.
}
}

\keyword{utilities}

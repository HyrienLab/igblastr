\name{auxdata-utils}

\alias{auxdata-utils}
\alias{auxdata_utils}
\alias{auxdata}

\alias{get_auxdata_path}
\alias{load_auxdata}
\alias{compute_germline_db_auxdata}

\alias{get_igblast_auxiliary_data}
\alias{load_igblast_auxiliary_data}

\title{Access IgBLAST auxiliary data}

\description{
  IgBLAST \emph{auxiliary data} is expected to annotate all the known
  germline J gene alleles for a given organism. It is provided by NCBI
  and is typically included in a standard IgBLAST installation.

  The \emph{auxiliary data} informs about the coding frame and CDR3/FWR4
  boundary on the J allele sequences. Although this information is not
  strictly needed by IgBLAST, the latter needs it in order to compute
  some of the AIRR fields like \code{vj_in_frame}, \code{productive},
  \code{cdr3}, \code{fwr4}, etc...

  \code{get_auxdata_path()} and \code{load_auxdata()} can be used
  to access the \emph{auxiliary data} included in IgBLAST or in
  one of the \emph{cached germline dbs} managed by \pkg{igblastr}.

  We also provide \code{compute_germline_db_auxdata()} to compute
  the \emph{auxiliary data} associated with the J alleles in a
  cached germline dbs managed by \pkg{igblastr}.
}

\usage{
## Access auxiliary data:
get_auxdata_path(organism, which=c("live", "original"))
load_auxdata(organism, which=c("live", "original"))

## Compute auxiliary data associated with the J alleles in germline db:
compute_germline_db_auxdata(db_name, ...)
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
  \item{db_name}{
    A single string specifying the name of a cached germline db.
    Use \code{\link{list_germline_dbs}()} to list all the cached
    germline dbs.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to
    \code{compute_auxdata()}. See \code{?\link{compute_auxdata}} for
    what arguments are supported.
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

  IMPORTANT NOTE: Using IgBLAST with \emph{auxiliary data} that is not
  compatible with the germline J gene sequences in the germline db can
  cause it to return improper frame status or CDR3 information (other
  returned information will still be correct).
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

  \code{compute_germline_db_auxdata()} returns the computed
  \emph{auxiliary data} in a data.frame with 1 row per J allele in
  the germline db and the same columns as the data.frame returned
  by \code{compute_auxdata()} or \code{read_auxdata()}.
  See \code{?\link{compute_auxdata}} for more information.
}

\seealso{
  \itemize{
    \item \code{\link{compute_auxdata}} to annotate a set of germline J
          gene allele sequences.

    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \code{\link{list_igblast_organisms}} to list the organisms
          supported by IgBLAST.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about IgBLAST \emph{auxiliary data}.

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
## Access IgBLAST auxiliary data for a given organism
## ---------------------------------------------------------------------

## IgBLAST only includes auxiliary data for the following organisms:
list_igblast_organisms()

get_auxdata_path("human")

human_auxdata <- load_auxdata("human")
head(human_auxdata)

## ---------------------------------------------------------------------
## Access the auxiliary data included in a germline db
## ---------------------------------------------------------------------

## List germline dbs with auxiliary data:
list_germline_dbs(with.auxdata.only=TRUE)

## Access the auxiliary data included in germline db
## _OGRDB.rhesus_monkey.IGH+IGK+IGL.202602 (this data
## was precomputed with compute_germline_db_auxdata()):
db_name <- "_OGRDB.rhesus_monkey.IGH+IGK+IGL.202602"

get_auxdata_path(db_name)

rhesus_monkey_auxdata <- load_auxdata(db_name)
head(rhesus_monkey_auxdata)

## ---------------------------------------------------------------------
## Compute the auxiliary data associated with the J alleles in a
## germline db
## ---------------------------------------------------------------------

## The auxiliary data included in germline db
## _OGRDB.rhesus_monkey.IGH+IGK+IGL.202602 was precomputed
## with compute_germline_db_auxdata(). Let's compute it again:
rhesus_monkey_auxdata2 <- compute_germline_db_auxdata(db_name)

## Sanity check:
stopifnot(identical(rhesus_monkey_auxdata2, rhesus_monkey_auxdata))

## Note that no auxiliary data was included in germline db
## _OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410 because it couldn't
## be properly computed for all the J alleles in the db:
db_name <- "_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410"
compute_germline_db_auxdata(db_name)  # warning!
}

\keyword{utilities}

\name{update_live_igdata}

\alias{update_live_igdata}
\alias{igdata_info}
\alias{time_since_live_igdata_last_checked}
\alias{reset_live_igdata}

\alias{print.auxdata_md5sum_df}
\alias{print.igdata_info}

\title{Update and manage IgBLAST auxiliary and internal data}

\description{
  A small set of low-level utility functions to update and manage
  IgBLAST \emph{auxiliary data} and \emph{internal data}.
}

\usage{
update_live_igdata(check.only=FALSE)

igdata_info()

time_since_live_igdata_last_checked(units="days")

reset_live_igdata(subdirs=c("all", "internal_data", "optional_file"))
}

\arguments{
  \item{check.only}{
    By default, \code{update_live_igdata()} checks for new IgBLAST
    auxiliary or internal data files available at NCBI, and it installs
    them if any are found. Set \code{check.only} to \code{TRUE} to only
    do the check without installing anything.
  }
  \item{units}{
    See \code{?base::\link{difftime}} for valid units.
  }
  \item{subdirs}{
    By default, \code{reset_live_igdata()} resets both \code{internal_data/}
    and \code{optional_file/} directories to their original states.
    Set \code{subdirs} to \code{"internal_data"} or \code{"optional_file"}
    to reset only a particular directory.
  }
}

\details{
  \subsection{Auxiliary data and internal data}{
    A standard IgBLAST installation -- like the one used by the
    \pkg{igblastr} package -- typically includes \emph{auxiliary data}
    and \emph{internal data} that are normally found in directories
    \code{internal_data/} and \code{optional_file/}, respectively.
    Both directories should be subdirectories of the \emph{root directory}
    of the IgBLAST installation, that is, of the directory returned by
    \code{\link{get_igblast_root}()}.

    We sometimes refer to this data simply as the \emph{IgBLAST data}.
  }

  \subsection{NCBI updates}{
    NCBI occasionally updates some of the files in the
    \code{internal_data/} and \code{optional_file/} directories between
    IgBLAST releases, and it is recommended to use the new files.
    They make the new files available at
    \url{https://ftp.ncbi.nih.gov/blast/executables/igblast/release/patch/}.

    To download and install these new files, simply call
    \code{update_live_igdata()}. This will check for new IgBLAST
    auxiliary or internal data files available at NCBI, and install
    them if any are found.

    You can restore the original files at any moment with
    \code{reset_live_igdata()}.
  }
}

\value{
  \code{update_live_igdata()} and \code{reset_live_igdata()} don't return
  anything (invisible \code{NULL}).

  \code{igdata_info()} returns a named list containing information about
  the state of the "live" and "original" IgBLAST data.

  \code{time_since_live_igdata_last_checked()} returns the time passed
  since the last run of \code{update_live_igdata()} in the specified
  units (days by default).
}

\seealso{
  \itemize{
    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

    \item \link{auxdata_utils} to access IgBLAST \emph{auxiliary data}.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \code{\link{get_igblast_root}} to get (or set) the IgBLAST
          installation used (or to be used) by the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## Check for NCBI updates
## ---------------------------------------------------------------------

igdata_info()
update_live_igdata(check.only=TRUE)
igdata_info()

## ---------------------------------------------------------------------
## "live" vs "original" IgBLAST data
## ---------------------------------------------------------------------

## By default, the "live IgBLAST data" gets accessed or returned:
get_auxdata_path("human")
live_human_auxdata <- load_auxdata("human")

## Access the original IgBLAST data:
get_auxdata_path("human", which="original")
orig_human_auxdata <- load_auxdata("human", which="original")

## "live" and "original" IgBLAST data can differ if the former was
## updated with update_live_igdata(). Otherwise, they'll be the same:
identical(live_human_auxdata, orig_human_auxdata)
}

\keyword{utilities}

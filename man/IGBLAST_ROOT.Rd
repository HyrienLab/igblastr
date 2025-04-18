\name{IGBLAST_ROOT}

\alias{IGBLAST_ROOT}

\title{Use an external IgBLAST installation}

\description{
  Select the \emph{external} IgBLAST installation to use in \pkg{igblastr}
  in a persistent manner.
}

\details{
  The \pkg{igblastr} package can use 2 types of IgBLAST installation:
  \enumerate{
    \item Internal (a.k.a. igblastr-managed): refers to an
          installation obtained with \code{\link{install_igblast}()}.

    \item External: refers to an installation that is not managed by
          the \pkg{igblastr} package. This is usually an installation
          that was manually performed by you or a system administrator on
          your machine. It can be a system-wide installation or a per-user
          installation.
  }
  To use an \emph{external} installation of IgBLAST in \pkg{igblastr},
  set environment variable IGBLAST_ROOT to the path of the installation.
  Note that this must be the path to the \emph{root directory} of the
  IgBLAST installation, that is, to a directory with a \code{bin}
  subdirectory in it that has the \code{igblastn} and \code{igblastp}
  \emph{standalone executables} (on Windows these executables are files
  named \code{igblastn.exe} and \code{igblastp.exe}, respectively).

  This can be done within your current R session with
  \code{Sys.setenv(IGBLAST_ROOT="path/to/igblast_root")}
  for testing. However, this won't be remembered across R sessions.

  To set IGBLAST_ROOT in a persistent manner, define it outside R.
  The exact way to do this is OS-dependent e.g. on Linux and Mac you
  can define it in your user's \code{.profile} by adding the following
  line to it:
  \preformatted{
      export IGBLAST_ROOT="path/to/igblast_root"
  }
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \code{\link{igblast_info}} to collect basic information about
          the IgBLAST installation used by the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\keyword{misc}

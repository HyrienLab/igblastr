\name{IGBLAST_ROOT}

\alias{IGBLAST_ROOT}

\title{Use an external IgBlast installation}

\description{
  Set environment variable IGBLAST_ROOT to use an \emph{external}
  installation of IgBlast.
}

\details{
  The \pkg{igblastr} package can use 2 types of IgBlast installation:
  \enumerate{
    \item Internal (a.k.a. igblastr-managed): refers to an
          installation obtained with \code{\link{install_igblast}()}.

    \item External: refers to an installation that is not managed by
          the \pkg{igblastr} package. This is usually an installation
          that was manually performed by you or a system administrator on
          your machine. It can be a system-wide installation or a per-user
          installation.
  }
  To use an \emph{external} installation of IgBlast in \pkg{igblastr},
  set environment variable IGBLAST_ROOT to the path of the installation.
  Note that this must be the path to the \emph{root directory} of the
  IgBlast installation, that is, to a directory with a \code{bin}
  subdirectory in it that has the \code{igblastn} and \code{igblastp}
  commands (\code{igblastn.exe} and \code{igblastp.exe} on Windows).

  This can be done within your current R session with
  \code{Sys.setenv(IGBLAST_ROOT="path/to/igblast_root")}
  for testing. However, this won't be remembered when you quit R!

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
    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBlast installation.

    \item \code{\link{igblast_info}} to check the IgBlast installation
          used by the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}


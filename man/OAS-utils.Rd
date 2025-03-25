\name{OAS-utils}

\alias{OAS-utils}
\alias{OAS_utils}

\alias{read_OAS_csv_metadata}
\alias{read_OAS_csv}
\alias{extract_sequences_from_paired_OAS_df}

\alias{list_paired_OAS_studies}
\alias{list_paired_OAS_units}
\alias{download_paired_OAS_units}

\alias{extract_metadata_from_OAS_units}
\alias{extract_sequences_from_paired_OAS_units}

\title{Download and manipulate OAS data}

\description{
  Some utility functions to query the Observed Antibody Space database,
  a.k.a. OAS, and to download and manipulate data from OAS.

  OAS's homepage: \url{https://opig.stats.ox.ac.uk/webapps/oas/}

  Note that OAS has two databases: the "Unpaired Sequences" database
  and the "Paired Sequences" database. Some of the utilities documented
  in this man page only work on data coming from the latter.
}

\usage{
## Read metadata/data from a single OAS unit file:
read_OAS_csv_metadata(file)
read_OAS_csv(file, skip=1, ...)
extract_sequences_from_paired_OAS_df(df, add.prefix=FALSE)

## Basic query of OAS website:
list_paired_OAS_studies(as.df=FALSE, recache=FALSE)
list_paired_OAS_units(study, as.df=FALSE, recache=FALSE)
download_paired_OAS_units(study, units=NULL, destdir=".", ...)

## Read metadata/data from a batch of downloaded OAS unit files:
extract_metadata_from_OAS_units(dir=".", pattern="\\\\.csv\\\\.gz$")
extract_sequences_from_paired_OAS_units(dir=".", pattern="\\\\.csv\\\\.gz$")
}

\arguments{
  \item{file}{
    A single string containing the path to an \emph{OAS unit file}.
  }
  \item{skip}{
    The number of lines of the data file to skip before beginning to read data.
    The first line in an OAS unit file contains metadata in JSON format,
    so must always be skipped.
  }
  \item{...}{
    For \code{read_OAS_csv()}: Extra arguments to be passed to the
    internal call to \code{read.table()}.
    See \code{?\link[utils]{read.table}} in the \pkg{utils} package for
    more information.

    For \code{download_paired_OAS_units()}: Extra arguments to be
    passed to the internal call to \code{download.file()}.
    See \code{?\link[utils]{download.file}} in the \pkg{utils} package for
    more information.
  }
  \item{df}{
    The data.frame or tibble returned by \code{read_OAS_csv()}.
  }
  \item{add.prefix}{
    \code{TRUE} or \code{FALSE}. Should the names on the returned
    \link[Biostrings]{DNAStringSet} object be the original sequence ids
    as-is (this is the default), or should the \code{heavy_chain_} and
    \code{light_chain_} prefixes be added to them?

    \code{extract_sequences_from_paired_OAS_df()} returns a
    \link[Biostrings]{DNAStringSet} object with the sequence ids as
    names. The sequence ids are obtained from the \code{sequence_id_heavy}
    and \code{sequence_id_light} columns of the supplied data.frame or
    tibble. By default, they are propagated as-is to the DNAStringSet
    object, which makes it difficult to recognize which chain (heavy or
    light) the antibody sequences are coming from. Setting \code{add.prefix}
    to \code{TRUE} will add the \code{heavy_chain_} or \code{light_chain_}
    prefix to the names on the DNAStringSet object, hence making it easy
    to identify which chain a given antibody sequence is coming from.
  }
  \item{as.df}{
    \code{TRUE} or \code{FALSE}. By default, i.e. when \code{as.df}
    is \code{FALSE}, \code{list_paired_OAS_studies()} and
    \code{list_paired_OAS_units()} return the list of studies or units
    in a character vector. Alternatively you can set \code{as.df} to
    \code{TRUE} to get the list in a 3-column data.frame that contains
    a directory index as displayed at
    \url{https://opig.stats.ox.ac.uk/webapps/ngsdb/paired/}
    or at
    \url{https://opig.stats.ox.ac.uk/webapps/ngsdb/paired/Jaffe_2022/csv/}.
  }
  \item{recache}{
    \code{TRUE} or \code{FALSE}. \code{list_paired_OAS_studies()} and
    \code{list_paired_OAS_units()} both cache the information retrieved
    from OAS website for the duration of the R session (note that this
    caching is done in memory so it does not persist across sessions).
    Set \code{recache} to \code{TRUE} to force a new retrieval (and
    recaching) of the results.
  }
  \item{study}{
    A single string containing the name of a study as returned by
    \code{list_paired_OAS_studies()}.
  }
  \item{units}{
    \code{NULL}, or a character vector that must be a subset of
    \code{list_paired_OAS_units(study)} in which case the download
    will be restricted to these units only.
  }
  \item{destdir}{
    A single string that is the path to the directory where the OAS unit
    files are to be downloaded.
  }
  \item{dir}{
    A single string that is the path to a directory containing OAS unit
    files. This will typically be the same as \code{destdir} above if the
    unit files were downloaded with \code{download_paired_OAS_units()}.
  }
  \item{pattern}{
    Regular expression passed to the internal call to \code{list.files()}
    to obtain the list of OAS unit files located in \code{dir}. No reason
    to change this unless you know what you are doing.
  }
}

\details{
  OAS delivers data in the form of \emph{OAS unit files}. These files are
  typically obtained by running the \code{bulk_download.sh} script that OAS
  generates based on one's search criteria. They are compressed CSV
  (comma-separated values) files with the \code{.csv.gz} extension.

  OAS unit files can vary a lot in size: from only a few KB to 25 MB or more.

  The first line in an OAS unit file contains metadata in JSON format (which
  means that these files cannot strictly be considered CSV files).

  The CSV data is MiAIRR-compliant (see The "MiAIRR format" paper in
  the References section below).
}

\value{
  \code{read_OAS_csv_metadata()} extracts the metadata from the specified
  OAS unit file and returns it in a named list.

  \code{read_OAS_csv()} extracts the data from the specified OAS unit file
  and returns it in a \link[tibble]{tibble}. The tibble has 1 row per
  antibody sequence if the data is unpaired (i.e. comes from the "Unpaired
  Sequences" database), or 1 row per sequence pair if the data is paired (i.e.
  comes from the "Paired Sequences" database).

  \code{extract_sequences_from_paired_OAS_df()} returns the sequence pairs
  in a named \link[Biostrings]{DNAStringSet} object where the names are
  the sequence ids. See \code{add.prefix} above for how the sequence ids
  are obtained.

  \code{list_paired_OAS_studies()} returns the list of studies that populate
  the "Paired Sequences" database in a character vector. This list
  can be seen here: \url{https://opig.stats.ox.ac.uk/webapps/ngsdb/paired/}.

  \code{list_paired_OAS_units()} returns the list of all the OAS unit files
  that belong to a given study from the "Paired Sequences" database.

  \code{download_paired_OAS_units()} returns an invisible \code{NULL}.

  \code{extract_metadata_from_OAS_units()} returns the metadata of all the
  OAS unit files found in the specified directory in a data.frame with 1 row
  per file.

  \code{extract_sequences_from_paired_OAS_units()} extracts the sequence
  pairs from all the OAS unit files found in the specified directory and
  returns them in a named \link[Biostrings]{DNAStringSet} object where the
  names are the sequence ids. The sequence ids are obtained by prefixing
  the original sequence ids found in the files with the name of the unit
  followed by \code{_heavy_chain_} or \code{_light_chain_}.
}

\references{
  \itemize{
    \item The OAS paper:

          Tobias H. Olsen, Fergus Boyles, Charlotte M. Deane.
          Observed Antibody Space: A diverse database of cleaned, annotated,
          and translated unpaired and paired antibody sequences.
          Protein Science (2021).
          \url{https://doi.org/10.1002/pro.4205}

    \item The "MiAIRR format" paper:

          Rubelt, F., Busse, C., Bukhari, S. et al. Adaptive Immune Receptor
          Repertoire Community recommendations for sharing immune-repertoire
          sequencing data. Nat Immunol 18, 1274–1278 (2017).
          \url{https://doi.org/10.1038/ni.3873}
  }
}

\seealso{
  \itemize{
    \item OAS's homepage at: \url{https://opig.stats.ox.ac.uk/webapps/oas/}

    \item The \code{\link{igblastn}} function to run the \emph{igblastn}
          executable from R. This is the main function in the \pkg{igblastr}
          package.

    \item \link[tibble]{tibble} objects implemented in the
          \pkg{tibble} package.

    \item \link[Biostrings]{DNAStringSet} objects implemented in the
          \pkg{Biostrings} package.
  }
}

\examples{
list_paired_OAS_studies()

list_paired_OAS_units("Eccles_2020")

## Import all the pairs of antibody sequences from the Eccles_2020 study:

download_dir <- tempdir()
download_paired_OAS_units("Eccles_2020", destdir=download_dir)

metadata <- extract_metadata_from_OAS_units(download_dir)
metadata  # data.frame with 1 row per unit file

sequences <- extract_sequences_from_paired_OAS_units(download_dir)
sequences  # DNAStringSet object

## Odd indices correspond to heavy chain sequences and even indices
## to light chain sequences:

head(names(sequences))

sequences[1:2]  # 1st pair
sequences[3:4]  # 2nd pair
sequences[5:6]  # 3rd pair
# etc...
}

\keyword{utilities}

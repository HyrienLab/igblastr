\name{download_OGRDB_germline_json}

\alias{download_OGRDB_germline_json}
\alias{extract_intdata_from_ogrdb_json}
\alias{extract_auxdata_from_ogrdb_json}

\title{Download AIRR-C JSON files from OGRDB}

\description{
  The \code{download_OGRDB_germline_json()} function downloads the
  AIRR-C JSON files assciated with the specified germline sets from
  OGRDB.

  The \code{extract_intdata_from_ogrdb_json()} function extracts the
  V gene delineations for the V alleles in an AIRR-C JSON file.
  Note that it's a simple reimplementation in R of the Python script
  \code{extract_intdata_from_ogrdb_json.py} included in IgBLAST.

  The \code{extract_auxdata_from_ogrdb_json()} function extracts the
  coding frame and CDR3 end information for the J alleles in an
  AIRR-C JSON file.

  If you use OGRDB in your research, please cite:
  \preformatted{  The current landscape of adaptive immune receptor genomic and
  repertoire data: OGRDB and VDJbase
  Lees, Peres, Klein, Amos et al., Nucleic Acids Research, November 2025.
  https://doi.org/10.1093/nar/gkaf1094}
}

\usage{
download_OGRDB_germline_json(organism, germline_sets,
                             source_set=FALSE,
                             destdir=".", overwrite=FALSE,
                             recache=FALSE, ...)

extract_intdata_from_ogrdb_json(json_path)
extract_auxdata_from_ogrdb_json(json_path)
}

\arguments{
  \item{organism}{
    A single string specifying the latin name of the organism for which
    to download the AIRR-C JSON file, e.g. \code{"Homo sapiens"}.
    See \url{https://ogrdb.airr-community.org/germline_sets} for the
    organisms that are currently available at OGRDB.
  }
  \item{germline_sets}{
    A named integer vector. The names on the vector must be the names of
    the germline sets to download and the values in the vector must be
    their versions.
    Go to \url{https://ogrdb.airr-community.org/germline_sets}, select a
    Species, and click on "Show Germline Sets" to see the list of germline
    sets available for that species and their versions.
  }
  \item{source_set}{
    Whether to download the AIRR-C JSON files assciated with the
    \emph{Source Sets} instead of the \emph{Reference Sets}. Only
    applies when \code{organism} is Homo sapiens.
  }
  \item{destdir}{
    A single string that is the path to the directory where the
    AIRR-C JSON files should be downloaded.
  }
  \item{overwrite}{
    \code{download_OGRDB_germline_json()} won't replace existing
    files, unless \code{overwrite} is set to \code{TRUE}.
  }
  \item{recache}{
    \code{download_OGRDB_germline_json()} uses a caching mechanism for
    the data that gets downloaded from OGRDB. The data is cached on disk so
    the caching is persistent across R sessions.

    Set \code{recache} to \code{TRUE} to force a new download (and recaching)
    of the data.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?\link[utils]{download.file}} in the \pkg{utils} package for
    more information.
  }
  \item{json_path}{
    The path to an AIRR-C JSON file as one obtained with
    \code{download_OGRDB_germline_json()}.
  }
}

\value{
  \code{download_OGRDB_germline_json()} returns the list of AIRR-C JSON
  files that it produced in an invisible character vector that carries
  the names of the corresponding germline sets.

  \code{extract_intdata_from_ogrdb_json()} returns the V gene delineations
  in a data.frame similar to the data.frame returned by \code{load_intdata()}.
  See \code{?\link{load_intdata}} for the details.

  \code{extract_auxdata_from_ogrdb_json()} returns the coding frame and
  CDR3 end information in a data.frame similar to the data.frame returned
  by \code{load_auxdata()}. See \code{?\link{load_auxdata}} for the details.
}

\seealso{
  \itemize{
    \item \code{\link{download_OGRDB_germline_sequences}} to download
          germline sequences from OGRDB.

    \item \code{\link{write_ndm_data}} to write V gene allele
          FWR/CDR annotations in \emph{ndm} format.

    \item \code{\link{compute_V_gene_delineations}} to annotate a set
          of germline V gene allele sequences.

    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

    \item The OGRDB website: \url{https://ogrdb.airr-community.org/}.
  }
}

\examples{
## ---------------------------------------------------------------------
## DOWNLOAD AIRR-C JSON FILES FOR MOUSE STRAIN PWD/PhJ FROM OGRDB
## ---------------------------------------------------------------------
## Visit https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
## for the list of germline sets available for mouse.

## Download AIRR-C JSON files for all BCR loci/regions for mouse strain
## PWD/PhJ to a temporary directory:
germline_sets <- c(`PWD/PhJ IGH`=2, `PWD/PhJ IGKV`=1, `PWD/PhJ IGLV`=1,
                   `IGKJ (all strains)`=1, `IGLJ (all strains)`=1)

dir.create(destdir <- tempfile("json_dir_"))
download_OGRDB_germline_json("Mus musculus", germline_sets,
                             destdir=destdir)

## List the downloaded files:
list.files(destdir)  # 5 AIRR-C JSON files (one per germline set)

## ---------------------------------------------------------------------
## EXTRACT V GENE DELINEATIONS FROM THE AIRR-C JSON FILES
## ---------------------------------------------------------------------
## It only makes sense to extract the V gene delineations from the
## AIRR-C JSON files that correspond to germline sets with V alleles,
## that is, from files IGH.json, IGKV.json, and IGLV.json:

json_path <- file.path(destdir, "IGH.json")
IGHV_intdata <- extract_intdata_from_ogrdb_json(json_path)
head(IGHV_intdata)

json_path <- file.path(destdir, "IGKV.json")
IGKV_intdata <- extract_intdata_from_ogrdb_json(json_path)

json_path <- file.path(destdir, "IGLV.json")
IGLV_intdata <- extract_intdata_from_ogrdb_json(json_path)

## Combine the 3 data.frames:
intdata <- rbind(IGHV_intdata, IGKV_intdata, IGLV_intdata)

## As a sanity check, we compare 'intdata' with the "internal data"
## included in germline db _OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410:
intdata0 <- load_intdata("_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410")

## We use same_ndm_data(), a small undocumented helper function, to
## verify that 'intdata' and 'intdata0' contain the same data (possibly
## with their rows in different order):
stopifnot(same_ndm_data(intdata, intdata0))

## Remove temporary directory:
unlink(destdir, recursive=TRUE)
}

\keyword{utilities}

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
  V gene delineations for the V alleles annotated in an AIRR-C JSON file.
  Note that it's a simple reimplementation in R of the Python script
  \code{makeogrannote.py} included in IgBLAST.

  The \code{extract_auxdata_from_ogrdb_json()} function extracts the
  coding frame and CDR3 end information for the J alleles annotated
  in an AIRR-C JSON file.

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

extract_intdata_from_ogrdb_json(json_path, extra_fields=NULL)
extract_auxdata_from_ogrdb_json(json_path, extra_fields=NULL)
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
  \item{extra_fields}{
    \code{NULL} (the default) or a character vector of valid AIRR-C
    "AlleleDescription fields" to extract from the AIRR-C JSON file
    and to include in the returned data.frame as additional columns. See
    \url{https://docs.airr-community.org/en/latest/datarep/germline.html#alleledescription-fields}
    for the list of valid "AlleleDescription fields".
  }
}

\value{
  \code{download_OGRDB_germline_json()} returns the list of AIRR-C JSON
  files that it produced in an invisible character vector that carries
  the names of the corresponding germline sets.

  \code{extract_intdata_from_ogrdb_json()} returns the V gene delineations
  in a data.frame with 1 row per V allele in the AIRR-C JSON file and
  with the same columns as the data.frame returned by \code{load_intdata()}.
  See \code{?\link{load_intdata}} for the details.
  If \code{extra_fields} was specified, the returned data.frame will
  have 1 additional column per field in \code{extra_fields}.

  \code{extract_auxdata_from_ogrdb_json()} returns the coding frame and
  CDR3 end information in a data.frame with 1 row per J allele in the
  AIRR-C JSON file and with the same columns as the data.frame returned by
  \code{load_auxdata()}. See \code{?\link{load_auxdata}} for the details.
  If \code{extra_fields} was specified, the returned data.frame will
  have 1 additional column per field in \code{extra_fields}.
}

\seealso{
  \itemize{
    \item \code{\link{download_OGRDB_germline_sequences}} to download
          germline sequences from OGRDB.

    \item \code{\link{write_ndm_data}} to write the data.frame returned
          by \code{extract_intdata_from_ogrdb_json()} in \emph{ndm} format.

    \item \code{\link{write_auxdata}} to write the data.frame returned
          by \code{extract_auxdata_from_ogrdb_json()} in \emph{.aux} format.

    \item \code{\link{compute_V_gene_delineations}} to annotate a set
          of germline V gene allele sequences.

    \item \code{\link{compute_auxdata}} to annotate a set
          of germline J gene allele sequences.

    \item \link{intdata_utils} to access IgBLAST \emph{internal data}.

    \item \link{auxdata_utils} to access IgBLAST \emph{auxiliary data}.

    \item The OGRDB website: \url{https://ogrdb.airr-community.org/}.
  }
}

\examples{
## ---------------------------------------------------------------------
## download_OGRDB_germline_json()
## ---------------------------------------------------------------------
## We're going to use download_OGRDB_germline_json() to download the
## AIRR-C JSON files for mouse strain PWD/PhJ from OGRDB.
## Visit https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
## for the list of germline sets available for mouse.

## We need the five following germline sets in order to cover all the
## BCR loci/regions for mouse strain PWD/PhJ:
germline_sets <- c(`PWD/PhJ IGH`=2, `PWD/PhJ IGKV`=1, `PWD/PhJ IGLV`=1,
                   `IGKJ (all strains)`=1, `IGLJ (all strains)`=1)

## Let's download them to a temporary directory:
dir.create(temp_json_dir <- tempfile("json_dir_"))
download_OGRDB_germline_json("Mus musculus", germline_sets,
                             destdir=temp_json_dir)

## List the downloaded files:
list.files(temp_json_dir)  # 5 AIRR-C JSON files (one per germline set)

## ---------------------------------------------------------------------
## extract_intdata_from_ogrdb_json()
## ---------------------------------------------------------------------
## Let's use extract_intdata_from_ogrdb_json() to extract the V gene
## delineations for the V alleles annotated in the AIRR-C JSON files.
## Note that germline sets "IGKJ (all strains)" and "IGLJ (all strains)"
## don't annotate any V alleles so we exclude files IGKJ.json and
## IGLJ.json:

json_path <- file.path(temp_json_dir, "IGH.json")
IGHV_intdata <- extract_intdata_from_ogrdb_json(json_path)
head(IGHV_intdata)

json_path <- file.path(temp_json_dir, "IGKV.json")
IGKV_intdata <- extract_intdata_from_ogrdb_json(json_path)
head(IGKV_intdata)

json_path <- file.path(temp_json_dir, "IGLV.json")
IGLV_intdata <- extract_intdata_from_ogrdb_json(json_path)
IGLV_intdata

## Combine the 3 data.frames:
intdata <- rbind(IGHV_intdata, IGKV_intdata, IGLV_intdata)
dim(intdata)

## 'intdata' should be the same as the "internal data" included in
## germline db _OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410:
intdata0 <- load_intdata("_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410")

## Let's use same_alleles_annot(), a small undocumented helper function,
## to verify that 'intdata' and 'intdata0' contain the same data (possibly
## with their rows in different order):
stopifnot(same_alleles_annot(intdata, intdata0))

## Note that calling extract_intdata_from_ogrdb_json() on an AIRR-C JSON
## file with no V alleles returns a 0-row data.frame with a warning:
json_path <- file.path(temp_json_dir, "IGKJ.json")
IGKJ_intdata <- extract_intdata_from_ogrdb_json(json_path)  # warning!
stopifnot(nrow(IGKJ_intdata) == 0)
json_path <- file.path(temp_json_dir, "IGLJ.json")
IGLJ_intdata <- extract_intdata_from_ogrdb_json(json_path)  # warning!
stopifnot(nrow(IGLJ_intdata) == 0)

## ---------------------------------------------------------------------
## extract_auxdata_from_ogrdb_json()
## ---------------------------------------------------------------------
## Let's use extract_auxdata_from_ogrdb_json() to extract the coding
## frame and CDR3 end information for the J alleles annotated in the
## AIRR-C JSON files. Note that germline sets "PWD/PhJ IGKV" and
## "PWD/PhJ IGLV" don't annotate any J alleles so we exclude files
## IGKV.json and IGLV.json:

json_path <- file.path(temp_json_dir, "IGH.json")
IGHJ_auxdata <- extract_auxdata_from_ogrdb_json(json_path)
IGHJ_auxdata

json_path <- file.path(temp_json_dir, "IGKJ.json")
IGKJ_auxdata <- extract_auxdata_from_ogrdb_json(json_path)
IGKJ_auxdata

json_path <- file.path(temp_json_dir, "IGLJ.json")
IGLJ_auxdata <- extract_auxdata_from_ogrdb_json(json_path)
IGLJ_auxdata

## Combine the 3 data.frames:
auxdata <- rbind(IGHJ_auxdata, IGKJ_auxdata, IGLJ_auxdata)
dim(auxdata)

## 'auxdata' should be the same as the "auxiliary data" included
## in germline db _OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410:
auxdata0 <- load_auxdata("_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410")
stopifnot(same_alleles_annot(auxdata, auxdata0))

## Note that calling extract_auxdata_from_ogrdb_json() on an AIRR-C JSON
## file with no J alleles returns a 0-row data.frame with a warning:
json_path <- file.path(temp_json_dir, "IGKV.json")
IGKV_auxdata <- extract_auxdata_from_ogrdb_json(json_path)  # warning!
stopifnot(nrow(IGKV_auxdata) == 0)
json_path <- file.path(temp_json_dir, "IGLV.json")
IGLV_auxdata <- extract_auxdata_from_ogrdb_json(json_path)  # warning!
stopifnot(nrow(IGLV_auxdata) == 0)

## Remove temporary directory:
unlink(temp_json_dir, recursive=TRUE)
}

\keyword{utilities}

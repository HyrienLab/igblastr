### Run this code **in this folder** to generate all the FASTA files in the
### various subfolders of this folder.

library(igblastr)

HUMAN_GERMLINE_SETS <- list(
    `202309`=c(IGH_VDJ=7, IGKappa_VJ=2, IGLambda_VJ=1),
    `202401`=c(IGH_VDJ=8, IGKappa_VJ=3, IGLambda_VJ=2),
    `202410`=c(IGH_VDJ=9, IGKappa_VJ=4, IGLambda_VJ=3)
)

download_human_germline_sequences <- function(overwrite=FALSE)
{
    expected_fasta_files <- paste0(igblastr:::IG_GROUPS, ".fasta")
    for (i in seq_along(HUMAN_GERMLINE_SETS)) {
        version <- names(HUMAN_GERMLINE_SETS)[[i]]
        germline_sets <- HUMAN_GERMLINE_SETS[[i]]
        for (flavor in c("ref", "src")) {
            source_set <- flavor == "src"
            destdir <- file.path(version, flavor)
            message("Creating FASTA files in ", destdir, " ... ",
                    appendLF=FALSE)
            filenames <- download_OGRDB_germline_sequences("Homo sapiens",
                                        germline_sets, source_set=source_set,
                                        destdir=destdir, overwrite=overwrite)
            stopifnot(length(filenames) == length(expected_fasta_files),
                      setequal(filenames, expected_fasta_files))
            message("ok")
        }
    }
}

make_human_auxdata <- function(version, germline_sets)
{
    dir.create(jsondir <- tempfile())
    for (flavor in c("ref", "src")) {
        source_set <- flavor == "src"
        destdir <- file.path(version, flavor)
        message("Creating IG[HKL]J_gl.aux files in ", destdir, " ... ",
                appendLF=FALSE)
        filenames <- download_OGRDB_germline_json("Homo sapiens",
                                    germline_sets, source_set=source_set,
                                    destdir=jsondir, overwrite=TRUE)
        stopifnot(identical(names(filenames), names(germline_sets)))
        for (filename in filenames) {
            json_path <- file.path(jsondir, filename)
            auxdata <- extract_auxdata_from_ogrdb_json(json_path)
            locus <- substr(filename, 1L, 3L)
            destfile <- file.path(destdir, paste0(locus, "J_gl.aux"))
            write_auxdata(auxdata, destfile)
        }
        message("ok")
    }
}

download_human_germline_sequences()

### Extract auxdata from the OGRDB json files but only for version 202410.
### Older versions break extract_auxdata_from_ogrdb_json()!
make_human_auxdata(names(HUMAN_GERMLINE_SETS)[[3]], HUMAN_GERMLINE_SETS[[3]])


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_human_intdata()
###
### Validates the intdata associated with the human germline sets listed
### in HUMAN_GERMLINE_SETS by comparing the two methods of acquisition:
###   1. Infer intdata from the gaps in the V allele sequences.
###   2. Extract intdata from OGRDB json file.
###
### WARNING: It looks like prior to version 3 (published in Oct 2024 and
### still the latest version as of March 2026), germline set IGLambda_VJ
### has inconsistent internal data!

validate_human_intdata <- function()
{
    organism <- "Homo sapiens"
    for (i in seq_along(HUMAN_GERMLINE_SETS)) {
        germline_sets <- HUMAN_GERMLINE_SETS[[i]]
        for (flavor in c("ref", "src")) {
            source_set <- flavor == "src"
            for (j in seq_along(germline_sets)) {
                germline_set <- germline_sets[j]
                what <- c(organism, " ", flavor, " germline set ",
                          "\"", names(germline_set), "\" ",
                          "(version ", germline_set, ")")
                message("Validating intdata for ", what, " ... ",
                        appendLF=FALSE)
                ok <- igblastr:::validate_OGRDB_intdata(organism, germline_set,
                                                        source_set=source_set)
                msg <- if (ok) "ok" else "DATA IS INCONSISTENT!"
                message(msg)
            }
        }
    }
}


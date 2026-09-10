### Run this code **in this folder** to generate all the FASTA files in the
### various subfolders of this folder.

library(igblastr)

RAINBOW_TROUT_GERMLINE_SETS <- list(
    `202606`=c(IGH_VDJ=1)
)

download_rainbow_trout_germline_sequences <- function(overwrite=FALSE)
{
    #expected_fasta_files <- paste0(igblastr:::IG_GROUPS, ".fasta")
    expected_fasta_files <- paste0("IGH", c("V", "D", "J"), ".fasta")
    for (i in seq_along(RAINBOW_TROUT_GERMLINE_SETS)) {
        destdir <- names(RAINBOW_TROUT_GERMLINE_SETS)[[i]]
        germline_sets <- RAINBOW_TROUT_GERMLINE_SETS[[i]]
        message("Creating FASTA files in ", destdir, " ... ",
                appendLF=FALSE)
        filenames <- download_OGRDB_germline_sequences("Oncorhynchus mykiss",
                                    germline_sets,
                                    destdir=destdir, overwrite=overwrite)
        stopifnot(length(filenames) == length(expected_fasta_files),
                  setequal(filenames, expected_fasta_files))
        message("ok")
    }
}

download_rainbow_trout_germline_sequences()

make_rainbow_trout_auxdata <- function(version, germline_sets, overwrite=FALSE)
{
    dir.create(json_dir <- tempfile())
    message("Creating IG[HKL]J_gl.aux files in ", version, " ... ",
            appendLF=FALSE)
    json_files <- download_OGRDB_germline_json("Oncorhynchus mykiss",
                                 germline_sets,
                                 destdir=json_dir, overwrite=TRUE)
    stopifnot(identical(names(json_files), names(germline_sets)))
    auxdata_files <-
        make_auxdata_files_from_ogrdb_jsons(json_dir, destdir=version,
                                            overwrite=overwrite)
    #expected_files <- sprintf("IG%sJ_gl.aux", c("H", "K", "L"))
    expected_files <- "IGHJ_gl.aux"
    stopifnot(identical(auxdata_files, expected_files))
    message("ok")
}

### Extract auxdata from the OGRDB json files.
for (i in seq_along(RAINBOW_TROUT_GERMLINE_SETS)) {
    version <- names(RAINBOW_TROUT_GERMLINE_SETS)[[i]]
    germline_sets <- RAINBOW_TROUT_GERMLINE_SETS[[i]]
    make_rainbow_trout_auxdata(version, germline_sets)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_rainbow_trout_intdata()
###
### Validates the intdata associated with the rainbow trout germline sets
### listed in RAINBOW_TROUT_GERMLINE_SETS by comparing the two methods of
### acquisition:
###   1. Infer intdata from the gaps in the V allele sequences.
###   2. Extract intdata from OGRDB json file.
###
### Note that all rainbow trout germline sets have consistent internal data.

validate_rainbow_trout_intdata <- function()
{
    organism <- "Oncorhynchus mykiss"
    for (i in seq_along(RAINBOW_TROUT_GERMLINE_SETS)) {
        germline_sets <- RAINBOW_TROUT_GERMLINE_SETS[[i]]
        for (j in seq_along(germline_sets)) {
            germline_set <- germline_sets[j]
            what <- c(organism, " germline set ",
                      "\"", names(germline_set), "\" ",
                      "(version ", germline_set, ")")
            message("Validating intdata for ", what, " ... ",
                    appendLF=FALSE)
            ok <- igblastr:::validate_OGRDB_intdata(organism, germline_set,
                      fwrcdr_widths=igblastr:::IMGT_RAINBOW_TROUT_FWRCDR_WIDTHS)
            msg <- if (ok) "ok" else "DATA IS INCONSISTENT!"
            message(msg)
        }
    }
}

validate_rainbow_trout_intdata()


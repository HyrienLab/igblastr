### Run this code **in this folder** to generate all the FASTA files in the
### strain-specific subfolders.

library(igblastr)

MOUSE_STRAIN_GERMLINE_SETS <- list(
    CAST_EiJ  =list(`202205`=c(IGH=1, IGKV=1, IGLV=1),
                    `202603`=c(IGH=2, IGKV=1, IGLV=2)),
    LEWES_EiJ =list(`202205`=c(IGH=1, IGKV=1, IGLV=1),
                    `202603`=c(IGH=2, IGKV=1, IGLV=1)),
    MSM_MsJ   =list(`202205`=c(IGH=1, IGKV=1, IGLV=1),
                    `202603`=c(IGH=1, IGKV=1, IGLV=2)),
    NOD_ShiLtJ=list(`202205`=c(IGH=1, IGKV=1, IGLV=1)),
    PWD_PhJ   =list(`202410`=c(IGH=2, IGKV=1, IGLV=1))
)

download_mouse_germline_sequences <- function(overwrite=FALSE)
{
    shared_germline_sets <- c(`IGKJ (all strains)`=1, `IGLJ (all strains)`=1)
    expected_files <- paste0(igblastr:::IG_GROUPS, ".fasta")
    for (i in seq_along(MOUSE_STRAIN_GERMLINE_SETS)) {
        strain_dir <- names(MOUSE_STRAIN_GERMLINE_SETS)[[i]]
        strain <- chartr("_", "/", strain_dir)
        strain_germline_sets <- MOUSE_STRAIN_GERMLINE_SETS[[i]]
        for (j in seq_along(strain_germline_sets)) {
            version <- names(strain_germline_sets)[[j]]
            message("Generating FASTA files for ", strain, " ",
                    "version ", version, " ... ", appendLF=FALSE)
            germline_sets <- strain_germline_sets[[j]]
            names(germline_sets) <- paste(strain, names(germline_sets))
            germline_sets <- c(germline_sets, shared_germline_sets)
            destdir <- file.path(strain_dir, version)
            filenames <- download_OGRDB_germline_sequences("Mus musculus",
                                        germline_sets,
                                        destdir=destdir, overwrite=overwrite)
            stopifnot(length(filenames) == length(expected_files),
                      setequal(filenames, expected_files))
            message("ok")
        }
    }
}

download_mouse_germline_sequences()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_mouse_intdata()
###
### Validates the intdata associated with the mouse germline sets listed
### in MOUSE_STRAIN_GERMLINE_SETS by comparing the two methods of
### acquisition:
###   1. Infer intdata from the gaps in the V allele sequences.
###   2. Extract intdata from OGRDB json file.
###
### WARNING: It looks like the following mouse germline sets have inconsistent
### internal data:
###
###   o "CAST/EiJ IGH"  version 1
###     Alleles with differences:
###     - IGHV0-NTLL*00: fwr3_end is 285 when inferred from the gaps,
###                      but is 283 in the JSON file;
###
###   o "CAST/EiJ IGLV" version 1
###     Alleles with differences:
###     - IGLV0-DUHW*00, IGLV0-EYCQ*00, IGLV0-JEYS*00, IGLV0-RQWY*00,
###       IGLV0-ZQBR*00: fwr3_end is 273 when inferred from the gaps,
###                      but is 245 in the JSON file;
###
###   o "LEWES/EiJ IGH" version 1
###     Alleles with differences:
###     - IGHV0-MJAE*00: fwr3_end is 288 when inferred from the gaps,
###                      but is 284 in the JSON file;
###
###   o "MSM/MsJ IGLV"  version 1
###     Alleles with differences:
###     - IGLV0-37NY*00, IGLV0-COOR*00: fwr3_end is 273 when inferred from
###                      the gaps, but is 245 in the JSON file.

validate_mouse_intdata <- function()
{
    organism <- "Mus musculus"
    for (i in seq_along(MOUSE_STRAIN_GERMLINE_SETS)) {
        strain_dir <- names(MOUSE_STRAIN_GERMLINE_SETS)[[i]]
        strain <- chartr("_", "/", strain_dir)
        strain_germline_sets <- MOUSE_STRAIN_GERMLINE_SETS[[i]]
        for (j in seq_along(strain_germline_sets)) {
            germline_sets <- strain_germline_sets[[j]]
            names(germline_sets) <- paste(strain, names(germline_sets))
            for (k in seq_along(germline_sets)) {
                germline_set <- germline_sets[k]
                what <- c(organism, " germline set ",
                          "\"", names(germline_set), "\" ",
                          "(version ", germline_set, ")")
                message("Validating intdata for ", what, " ... ",
                        appendLF=FALSE)
                ok <- igblastr:::validate_OGRDB_intdata(organism, germline_set)
                msg <- if (ok) "ok" else "DATA IS INCONSISTENT!"
                message(msg)
            }
        }
    }
}


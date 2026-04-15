### =========================================================================
### Parse and clean IMGT FASTA headers
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### parse_imgt_fasta_headers()
###

### The IMGT FASTA headers contain 15 fields separated by '|'.
### See https://www.imgt.org/IMGTindex/Fasta.php
.IMGT_FASTA_HEADERS <- c(
    IMGT_acc=
        "1. IMGT/LIGM-DB accession number(s)",
    allele_name=
        "2. IMGT gene and allele name",
    organism=
        "3. species",
    func=
        "4. IMGT allele functionality",
    region=
        "5. exon(s), region name(s), or extracted label(s)",
    startend_in_IMGT_acc=
        "6. start and end positions in the IMGT/LIGM-DB accession number(s)",
    nb_nuc=
        "7. number of nucleotides in the IMGT/LIGM-DB accession number(s)",
    codon_start=
        "8. codon start, or 'NR' (not relevant) for non coding labels",
    extra_nuc_5prime=
        paste0("9. +n: number of nucleotides (nt) added in 5' compared ",
               "to the corresponding label extracted from IMGT/LIGM-DB"),
    extra_nuc_3prime=
        paste0("10. +n or -n: number of nucleotides (nt) added or removed ",
               "in 3' compared to the corresponding label extracted ",
               "from IMGT/LIGM-DB"),
    nuc_corrected=
        paste0("11. +n, -n, and/or nS: number of added, deleted, and/or ",
               "substituted nucleotides to correct sequencing errors, ",
               "or 'not corrected' if non corrected sequencing errors"),
    nb_aa=
        paste0("12. number of amino acids (AA): this field indicates that ",
               "the sequence is in amino acids"),
    nb_chars=
        "13. number of characters in the sequence: nt (or AA)+IMGT gaps=total",
    partial=
        "14. partial (if it is)",
    revcomp=
        "15. reverse complementary (if it is)"
)

### Returns a 15-column character matrix with 1 row per header.
parse_imgt_fasta_headers <- function(headers)
{
    stopifnot(is.character(headers))
    if (anyNA(headers))
        stop(wmsg("some headers are NAs"))
    if (any(is_white_str(headers)))
        stop(wmsg("some headers are empty"))
    header_parts <- strsplit(headers, "|", fixed=TRUE)
    data <- right_pad_with_empty_strings_and_unlist(header_parts, width=15L)
    matrix(data, ncol=15L, byrow=TRUE,
           dimnames=list(names(headers), names(.IMGT_FASTA_HEADERS)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_imgt_fasta_headers()
###

### Not exported!
### Performs the same FASTA header cleanup as the edit_imgt_file.pl script
### included in IgBLAST. Note that the latter does some funky business with
### IG allele names for Mus spretus, so we mimick it here.
### Returns the germline gene allele names.
clean_imgt_fasta_headers  <- function(headers, what="some allele names")
{
    stopifnot(is.character(headers))
    if (anyNA(headers))
        stop(wmsg(what, " are NAs"))
    if (any(is_white_str(headers)))
        stop(wmsg(what, " are empty"))

    header_parts <- CharacterList(strsplit(headers, "|", fixed=TRUE))

    ## Extract 2nd field. Do NOT use parse_imgt_fasta_headers() for this.
    ## Note that this extraction method will get the 1st field if a header
    ## has no pipe or if it has only one pipe with nothing after it.
    ## Because of this, clean_imgt_fasta_headers() is idempotent i.e. it
    ## will be a no-op when called on headers that already went thru it.
    allele_names <- tails(heads(header_parts, n=2L), n=1L)
    stopifnot(all(lengths(allele_names) == 1L))
    allele_names <- trimws2(as.character(allele_names))
    if (!all(nchar(allele_names) >= 2L))
        stop(wmsg(what, " are less than 2-character long"))

    ## Implement funky business with Mus spretus: append _Mus_spretus suffix
    ## to IG allele names if species reported in 3rd field is Mus spretus.
    funky_idx <- which(substr(allele_names, 1L, 2L) == "IG" &
                       lengths(header_parts) >= 3L)
    if (length(funky_idx) != 0L) {
        ## Extract 3rd field.
        species <- tails(heads(header_parts[funky_idx], n=3L), n=1L)
        stopifnot(all(lengths(species) == 1L))
        idx <- grep("Mus\\s+spretus", as.character(species))
        if (length(idx) != 0L) {
            funky_idx <- funky_idx[idx]
            allele_names[funky_idx] <-
                paste0(allele_names[funky_idx], "_Mus_spretus")
        }
    }

    allele_names
}


\name{parse_imgt_fasta_headers}

\alias{parse_imgt_fasta_headers}

\title{Parse IMGT FASTA headers}

\description{
  The IMGT FASTA headers contain 15 fields separated by \code{|}.
  See \url{https://www.imgt.org/IMGTindex/Fasta.php}.

  \code{parse_imgt_fasta_headers()} takes a vector of headers and parse
  them into a 15-column matrix with 1 row per header.
}

\usage{
parse_imgt_fasta_headers(headers)
}

\arguments{
  \item{headers}{
    A character vector of IMGT FASTA headers.
  }
}

\value{
  A 15-column character matrix with 1 row per header. The column names on
  the matrix are:
  \enumerate{
    \item \code{IMGT_acc}: IMGT/LIGM-DB accession number(s);
    \item \code{allele_name}: IMGT gene and allele name;
    \item \code{organism}: species/organism;
    \item \code{func}: IMGT allele functionality;
    \item \code{region}: exon(s), region name(s), or extracted label(s);
    \item \code{startend_in_IMGT_acc}: start and end positions in the
          IMGT/LIGM-DB accession number(s);
    \item \code{nb_nuc}: number of nucleotides in the IMGT/LIGM-DB
          accession number(s);
    \item \code{codon_start}: codon start, or 'NR' (not relevant) for
          non coding labels;
    \item \code{extra_nuc_5prime}: +n: number of nucleotides (nt) added in 5'
          compared to the corresponding label extracted from IMGT/LIGM-DB;
    \item \code{extra_nuc_3prime}: +n or -n: number of nucleotides (nt) added
          or removed in 3' compared to the corresponding label extracted
          from IMGT/LIGM-DB;
    \item \code{nuc_corrected}: +n, -n, and/or nS: number of added, deleted,
          and/or substituted nucleotides to correct sequencing errors,
          or \code{not corrected} if non corrected sequencing errors;
    \item \code{nb_aa}: number of amino acids (AA): this field indicates
          that the sequence is in amino acids;
    \item \code{nb_chars}: number of characters in the sequence: nt
          (or AA)+IMGT gaps=total;
    \item \code{partial}: partial (if it is);
    \item \code{revcomp}: reverse complementary (if it is).
  }
}

\seealso{
  \itemize{
    \item \url{https://www.imgt.org/IMGTindex/Fasta.php} for the format
          of IMGT FASTA headers.

    \item \code{\link{download_IMGT_germline_sequences}} to download
          germline sequences from IMGT.

    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.
  }
}

\examples{
if (IMGT_is_up()) {
  ## -------------------------------------------------------------------
  ## Retrieve IMGT FASTA files for mouse and parse their headers
  ## -------------------------------------------------------------------

  ## As of April 10, 2026, the latest IMGT/V-QUEST release is 202614-2:
  list_IMGT_releases()

  list_IMGT_organisms("202614-2")

  ## Download Mouse BCR germline gene allele sequences from IMGT/V-QUEST
  ## 202614-2 to a temporary directory:
  dir.create(destdir <- tempfile("fasta_dir_"))
  download_IMGT_germline_sequences("202614-2", organism="Mus musculus",
                                   destdir=destdir)

  ## List the downloaded files:
  list.files(destdir)  # 7 FASTA files

  ## Load the V allele sequences (gapped):
  IMGT_groups <- c("IGHV", "IGKV", "IGLV")
  V_fasta_files <- file.path(destdir, paste0(IMGT_groups, ".fasta"))
  gapped_V_alleles <- readDNAStringSet(V_fasta_files)

  ## 'gapped_V_alleles' is a DNAStringSet object that carries the IMGT
  ## FASTA headers as names:
  gapped_V_alleles              # 865 V alleles
  names(gapped_V_alleles)[1:4]  # first 4 headers

  ## Let's parse the headers:
  parsed_V_headers <- parse_imgt_fasta_headers(names(gapped_V_alleles))
  dim(parsed_V_headers)    # 865 x 15 matrix
  parsed_V_headers[1:4, ]  # first 4 rows

  ## -------------------------------------------------------------------
  ## Scrutinize IMGT V alleles reported as "partial in 5'"
  ## -------------------------------------------------------------------

  ## The presence of "partial in 5'" in the 14th field ('partial') of
  ## a header indicates a sequence that is truncated at the 5' end:
  partial <- parsed_V_headers[ , "partial"]
  table(partial)
  is_partial_in_5prime <- grepl("in 5", partial)
  table(is_partial_in_5prime)  # 35 sequences marked as "partial in 5'"

  ## The corresponding sequences are expected to start with a gap:
  start_with_gap <- grepl("^\\\\.", as.character(gapped_V_alleles))

  ## However, it seems that 1 V allele (IGHV8-9*02) is marked
  ## as "partial in 5'" but does NOT start with a gap:
  table(is_partial_in_5prime, start_with_gap)
  not_ok <- is_partial_in_5prime & !start_with_gap
  parsed_V_headers[not_ok, ]  # partial in 5'
  gapped_V_alleles[not_ok]  # does NOT start with a gap!

  not_ok_V_allele <- unname(parsed_V_headers[not_ok, "allele_name"])

  stopifnot(identical(not_ok_V_allele, "IGHV8-9*02"))  # sanity check

  ## -------------------------------------------------------------------
  ## Scrutinize IMGT V alleles reported as "partial in 3'"
  ## -------------------------------------------------------------------

  ## The presence of "partial in 3'" in the 14th field ('partial') of
  ## a header indicates a sequence that is truncated at the 3' end:
  is_partial_in_3prime <- grepl("in 3", partial)
  table(is_partial_in_3prime)  # 72 sequences marked as "partial in 3'"

  ## Note that the FWR3 always ends at position 312 on a gapped V allele
  ## sequence. This means that any gapped sequence shorter than that should
  ## be considered truncated at the 3' end:
  gapped_seq_lens <- lengths(gapped_V_alleles)
  has_incomplete_fwr3 <- gapped_seq_lens < 312

  ## However, it seems that 12 V alleles have an incomplete FWR3 but
  ## are NOT marked as "partial in 3'":
  table(has_incomplete_fwr3, is_partial_in_3prime)
  not_ok <- has_incomplete_fwr3 & !is_partial_in_3prime
  gapped_V_alleles[not_ok]  # incomplete FWR3
  parsed_V_headers[not_ok, ]  # NOT marked as partial in 3'

  parsed_V_headers[not_ok, "allele_name"]  # suspect alleles

  stopifnot(sum(not_ok) == 12L)  # sanity check

  ## -------------------------------------------------------------------
  ## A quick look at the IMGT header of J alleles
  ## -------------------------------------------------------------------

  ## Load the J allele sequences:
  IMGT_groups <- c("IGHJ", "IGKJ", "IGLJ")
  J_fasta_files <- file.path(destdir, paste0(IMGT_groups, ".fasta"))
  J_alleles <- readDNAStringSet(J_fasta_files)

  ## 'J_alleles' is a DNAStringSet object that carries the IMGT
  ## FASTA headers as names:
  J_alleles              # 27 J alleles
  names(J_alleles)[1:4]  # first 4 headers

  ## Let's parse the headers:
  parsed_J_headers <- parse_imgt_fasta_headers(names(J_alleles))
  dim(parsed_J_headers)    # 27 x 15 matrix
  parsed_J_headers[1:4, ]  # first 4 rows

  ## Codon start is the start position of the first codon:
  as.integer(parsed_J_headers[ , "codon_start"])

  ## Remove temporary directory:
  unlink(destdir, recursive=TRUE)
}
}

\keyword{utilities}

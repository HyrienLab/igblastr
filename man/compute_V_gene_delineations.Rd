\name{compute_V_gene_delineations}

\alias{compute_V_gene_delineations}
\alias{compute_imgt_intdata}

\title{Compute V gene delineations}

\description{
  A utility function to compute the FWR/CDR boundaries (a.k.a. V gene
  delineations) for a set of germline V gene allele sequences.
  Note that the input sequences must be \emph{gapped}, that is, they
  must contain the so-called "IMGT gaps" whose purpose is to convey
  the delineation information.
}

\usage{
compute_V_gene_delineations(gapped_V_alleles, as.IRangesList=FALSE)
}

\arguments{
  \item{gapped_V_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline V
    gene allele \emph{gapped} sequences.
  }
  \item{as.IRangesList}{
    Experts only! Set to \code{TRUE} to return the delineations in an
    \link[IRanges]{IRangesList} object instead of a data.frame. In
    this case, the returned IRangesList object will have 1 list element
    per supplied V allele sequence and will carry the following metadata
    columns (accessible with \code{\link[S4Vectors]{mcols}()}):
    \code{seq_len}, \code{coding_frame_start}, \code{starting_gap},
    \code{all_gaps_in_frame}, \code{all_gaps_contained}.
    See Value section below for what these columns contain.
  }
}

\details{
  See
  \url{https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html}
  and
  \url{https://www.imgt.org/IMGTScientificChart/Numbering/IMGT-Kabat_part1.html}
  for how gaps are used in the context of "IMGT unique numbering for
  V-DOMAIN and V-LIKE-DOMAIN".

  Note that the annotation returned by \code{compute_V_gene_delineations()}
  can be used by IgBLAST as a substitute to the \emph{internal data} shipped
  with the IgBLAST software (and typically included in a standard IgBLAST
  installation).

  See \code{?\link{load_intdata}} for more information about IgBLAST
  \emph{internal data}.
}

\value{
  Returns the computed \emph{internal data} in a data.frame with 1 row
  per supplied V allele sequence and a set of columns similar to the
  data.frame returned by \code{\link{read_ndm_data}()}. With the following
  differences:
  \itemize{
    \item The data.frame returned by \code{compute_V_gene_delineations()} has
          no \code{chain_type} column.
    \item The data.frame returned by \code{compute_V_gene_delineations()} has
          the following additional columns:
          \enumerate{
            \item \code{seq_len}: length of the sequence after removal of
                  the gaps;
            \item \code{starting_gap}: size (in number of nucleotides) of
                  gap block located at the very beginning of the sequence
                  if any (0 if no such block);
            \item \code{all_gaps_in_frame}: indicates whether the gap blocks
                  align with the underlying coding frame or not;
            \item \code{all_gaps_contained}: \code{TRUE} if the gap blocks
                  don't cross the FWR/CDR boundaries, and \code{FALSE}
                  otherwise.
          }
  }
}

\seealso{
  \itemize{
    \item \url{https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html}
          for how the "IMGT unique numbering for V-DOMAIN and V-LIKE-DOMAIN"
          works.

    \item \code{\link{load_intdata}} to access IgBLAST \emph{internal data}.

    \item \code{\link{compute_auxdata}} to annotate a set of germline J
          gene allele sequences.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. USING HUMAN GERMLINE V GENE ALLELE SEQUENCES FROM OGRDB
## ---------------------------------------------------------------------
## Visit https://ogrdb.airr-community.org/germline_sets/Homo%20sapiens
## for the list of human germline sets available at OGRDB.

## Download OGRDB germline gene allele sequences for human locus IGK
## to a temporary directory:
germline_set <- c(`IGKappa_VJ`=4)

dir.create(tmpdir <- tempfile())
download_OGRDB_germline_sequences("Homo sapiens", germline_set,
                                  destdir=tmpdir)

## List the downloaded files:
list.files(tmpdir)  # 2 FASTA files (one per IMGT group)

## The V allele sequences are **gapped**:
IGKV_alleles <- readDNAStringSet(file.path(tmpdir, "IGKV.fasta"))
as.character(IGKV_alleles[1:3])

## Compute the V gene delineations for 'IGKV_alleles':
IGKV_gene_delinea <- compute_V_gene_delineations(IGKV_alleles)
head(IGKV_gene_delinea)

## As a sanity check, we're going to compare 'IGKV_gene_delinea' with
## the V gene delineations provided by the corresponding AIRR-C JSON
## file from OGRDB.

download_OGRDB_germline_json("Homo sapiens", germline_set,
                             destdir=tmpdir)

list.files(tmpdir)  # 1 JSON file

json_path <- file.path(tmpdir, "IGKappa_VJ.json")
IGKV_gene_delinea2 <- extract_intdata_from_ogrdb_json(json_path)
head(IGKV_gene_delinea2)

## Keep only the shared columns:
shared_colnames <- intersect(colnames(IGKV_gene_delinea),
                             colnames(IGKV_gene_delinea2))
IGKV_gene_delinea  <- IGKV_gene_delinea[ , shared_colnames]
IGKV_gene_delinea2 <- IGKV_gene_delinea2[ , shared_colnames]

## We use same_ndm_data(), a small undocumented helper function, to
## verify that 'IGKV_gene_delinea' and 'IGKV_gene_delinea2' contain
## the same delineation data (possibly with their rows in different
## order):
stopifnot(same_ndm_data(IGKV_gene_delinea, IGKV_gene_delinea2))

## Remove temporary directory:
unlink(tmpdir, recursive=TRUE)

## ---------------------------------------------------------------------
## 2. USING MOUSE GERMLINE V GENE ALLELE SEQUENCES FROM OGRDB
## ---------------------------------------------------------------------
## Visit https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
## for the list of mouse germline sets available at OGRDB.

## Download OGRDB germline V gene allele sequences for locus IGL of
## mouse strain MSM/MsJ to a temporary directory:
germline_set <- c(`MSM/MsJ IGLV`=1)

dir.create(tmpdir <- tempfile())
download_OGRDB_germline_sequences("Mus musculus", germline_set,
                                  destdir=tmpdir)

## List the downloaded files:
list.files(tmpdir)  # 1 FASTA file

## Compute the V gene delineations for 'IGKV_alleles':
IGLV_alleles <- file.path(tmpdir, "IGLV.fasta")
IGLV_gene_delinea <- compute_V_gene_delineations(IGLV_alleles)
IGLV_gene_delinea

## Note that for the first two alleles, the length of the **ungapped**
## sequence (seq_len) is smaller than the FWR3 end position (fwr3_end):
IGLV_gene_delinea[ , c("fwr3_end", "seq_len")]

## In other words, the sequences of alleles IGLV0-37NY*00 and
## IGLV0-COOR*00 end before the end of their FWR3.

## This is reflected in the V gene delineations provided by the
## corresponding AIRR-C JSON file from OGRDB:
download_OGRDB_germline_json("Mus musculus", germline_set,
                             destdir=tmpdir)
list.files(tmpdir)
json_path <- file.path(tmpdir, "IGLV.json")
IGLV_gene_delinea2 <- extract_intdata_from_ogrdb_json(json_path)
IGLV_gene_delinea2

## As you can see, when a V allele sequence ends before the end of
## its FWR3, the FWR3 end position reported in the AIRR-C JSON file
## is the position of the last nucleotide in the (ungapped) sequence.

## Remove temporary directory:
unlink(tmpdir, recursive=TRUE)
}

\keyword{utilities}
